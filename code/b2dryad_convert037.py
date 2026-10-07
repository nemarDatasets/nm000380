"""IEEG037 (Dryad doi:10.7272/Q6VD6WM2; Rao, Sellers et al. 2018 Curr Biol) -> iEEG-BIDS (stimulation + SPS recordings).

Converted: iEEG_DuringStimulation_EC*.zip (raw iEEG before/during/after continuous 100 Hz stimulation; Natus 'clinical',
TDT and NeuroOmega files) -> task-stim; iEEG_SPS.zip (single-pulse stimulation) -> task-sps.
Representation: BrainVision IEEE_FLOAT_32 multiplexed; the source float64 samples are rounded to float32 (max abs
rounding error per file recorded in code/conversion_report.json; the float64 originals are in sourcedata). No filtering,
resampling, re-referencing or channel removal.
NOT converted here: iEEG_NaturalBehavior.zip (ECoG segments around mood reports): its ECoG.time fields hold absolute
MATLAB datenums (real recording dates/times), so neither the zip nor its timestamps are published (kept private).
Usage: python b2dryad_convert037.py <src_dir> <bids_root>
"""
import csv, hashlib, json, os, re, shutil, sys, zipfile
import numpy as np, h5py, scipy.io as sio
try:
    import zipfile_deflate64  # noqa: F401
except ImportError:
    pass

SRC, OUT = sys.argv[1], sys.argv[2]
TMP = "/work/x037"
os.makedirs(OUT, exist_ok=True)
os.makedirs(TMP, exist_ok=True)


def wtsv(p, h, rows):
    with open(p, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(h)
        for r in rows:
            w.writerow(["n/a" if (v is None or v == "") else v for v in r])


def wjson(p, o):
    with open(p, "w") as f:
        json.dump(o, f, indent=2, ensure_ascii=False)
        f.write("\n")


def h5cell(f, ds):
    a = np.array(ds)
    out = np.empty(a.shape, dtype=object)
    for idx, r in np.ndenumerate(a):
        o = f[r]
        v = np.array(o).ravel()
        out[idx] = "".join(chr(int(c)) for c in v) if o.attrs.get("MATLAB_class", b"") == b"char" else (v.tolist() if v.size else "")
    return out


def elecs_v5(path):
    m = sio.loadmat(path, squeeze_me=False, struct_as_record=False)
    an = m.get("anatomy")
    em = m.get("elecmatrix")
    rows = []
    if an is not None and em is not None and an.shape[0] == em.shape[0]:
        for i in range(an.shape[0]):
            g = lambda j: str(np.squeeze(an[i, j])) if an.shape[1] > j else ""
            rows.append({"name": g(0), "long": g(1), "type": g(2), "region": g(3), "xyz": em[i].tolist()})
    return rows


def parse_name(fn):
    b = os.path.basename(fn)[:-4]
    amp = re.search(r"(\d+)mA", b)
    site = None
    for tok, nm in [("MedOFC", "medial OFC"), ("Insula", "insula"), ("VentralCingulate", "ventral cingulate"), ("DorsalCingulate", "dorsal cingulate"),
                    ("Area25", "area 25 (subgenual cingulate)"), ("SuperiorCingulate", "superior cingulate")]:
        if tok in b:
            site = nm
    side = "right" if "_Right" in b else ("left" if "_Left" in b else None)
    system = "neuroomega" if "_NO_" in b else ("tdt" if re.search(r"_B\d+_rawData|_\d+mA_rawData", b) else "natus")
    if site is None and system == "natus":
        site = "lateral OFC (README: no location in the file name = lateral OFC)"
    return {"amplitude_mA": int(amp.group(1)) if amp else None, "site": site, "side": side, "system": system, "stem": b}


report = {"runs": [], "skipped": []}
runs = {}
subjects = {}
zips = sorted(z for z in os.listdir(SRC) if z.startswith("iEEG_DuringStimulation_") and z.endswith(".zip")) + ["iEEG_SPS.zip"]
for zn in zips:
    Z = zipfile.ZipFile(os.path.join(SRC, zn))
    mats = [i for i in Z.infolist() if i.filename.endswith(".mat")]
    # electrode file of this subject (v5 clinical_elecs_all.mat or TDT_elecs_all.mat)
    el_rows = []
    for i in mats:
        if os.path.basename(i.filename) in ("clinical_elecs_all.mat", "TDT_elecs_all.mat"):
            p = Z.extract(i, TMP)
            try:
                el_rows = elecs_v5(p)
            except Exception as e:
                report["skipped"].append({"file": i.filename, "reason": "electrode file unreadable: " + repr(e)[:120]})
            os.remove(p)
    for i in mats:
        bn = os.path.basename(i.filename)
        if bn in ("clinical_elecs_all.mat", "TDT_elecs_all.mat", "TDT_elecs_special.mat", "SPS_Group.mat"):
            continue
        task = "sps" if zn == "iEEG_SPS.zip" else "stim"
        subj = re.match(r"(EC\d+)", bn).group(1)
        meta = parse_name(bn)
        p = Z.extract(i, TMP)
        h = hashlib.sha256()
        with open(p, "rb") as fh:
            for blk in iter(lambda: fh.read(1 << 24), b""):
                h.update(blk)
        src_sha = h.hexdigest()
        try:
            f = h5py.File(p, "r")
        except OSError as e:
            report["skipped"].append({"file": f"{zn}::{i.filename}", "reason": f"not readable as MATLAB v7.3/HDF5 ({e!r})"[:200], "sha256_extracted": src_sha})
            os.remove(p)
            print("SKIP", bn, flush=True)
            continue
        dkey = next((k for k in ("rawData", "RawData") if k in f), None) or next((k for k in f.keys() if k.startswith("during")), None)
        fs = float(np.array(f["Fs"]).ravel()[0])
        ds = f[dkey]
        ns, nch = ds.shape
        akey = next((k for k in ("selectAnatomy", "selectAnatmy", "Anatomy", "Anatmy") if k in f), None)
        names, longn, etype, region = [f"ch{k + 1:03d}" for k in range(nch)], [""] * nch, [""] * nch, [""] * nch
        names_src = "none in file (generic ch###)"
        if akey:
            c = h5cell(f, f[akey])
            if c.shape[1] == nch:
                names = [str(x) for x in c[0]]
                names_src = akey
                if c.shape[0] >= 4:
                    longn, etype, region = [str(x) for x in c[1]], [str(x) for x in c[2]], [str(x) for x in c[3]]
        # unique names (BIDS); keep originals in description
        seen = {}
        bnames = []
        for nm in names:
            nm2 = nm if nm not in ("", "NaN") else "unnamed"
            seen[nm2] = seen.get(nm2, 0) + 1
            bnames.append(nm2 if seen[nm2] == 1 else f"{nm2}_{seen[nm2]}")
        sub = f"sub-{subj}"
        subjects.setdefault(sub, {"el": el_rows, "systems": set()})
        if el_rows and not subjects[sub]["el"]:
            subjects[sub]["el"] = el_rows
        subjects[sub]["systems"].add(meta["system"])
        acq = meta["system"]
        k = runs.get((sub, task, acq), 0) + 1
        runs[(sub, task, acq)] = k
        stem = f"{sub}_task-{task}_acq-{acq}_run-{k:02d}"
        d = os.path.join(OUT, sub, "ieeg")
        os.makedirs(d, exist_ok=True)
        maxerr, absmax = 0.0, 0.0
        step = max(1, int(2e8 // (8 * nch)))
        with open(os.path.join(d, stem + "_ieeg.eeg"), "wb") as fo:
            for s0 in range(0, ns, step):
                x = ds[s0:s0 + step, :]
                x32 = x.astype("<f4")
                fin = np.isfinite(x)
                if fin.any():
                    maxerr = max(maxerr, float(np.max(np.abs(x32[fin].astype(np.float64) - x[fin]))))
                    absmax = max(absmax, float(np.max(np.abs(x[fin]))))
                fo.write(x32.tobytes())
        unit = "V" if meta["system"] == "tdt" else "µV"
        vh = ["Brain Vision Data Exchange Header File Version 1.0", f"; Written by b2dryad_convert037.py from {zn}::{i.filename} (float64 rounded to float32)", "",
              "[Common Infos]", "Codepage=UTF-8", f"DataFile={stem}_ieeg.eeg", f"MarkerFile={stem}_ieeg.vmrk", "DataFormat=BINARY", "DataOrientation=MULTIPLEXED",
              f"NumberOfChannels={nch}", f"SamplingInterval={1e6 / fs!r}", "", "[Binary Infos]", "BinaryFormat=IEEE_FLOAT_32", "", "[Channel Infos]"]
        vh += [f"Ch{j + 1}={nm.replace(',', chr(92) + '1')},,1,{unit}" for j, nm in enumerate(bnames)]
        open(os.path.join(d, stem + "_ieeg.vhdr"), "w", encoding="utf-8").write("\n".join(vh) + "\n")
        open(os.path.join(d, stem + "_ieeg.vmrk"), "w", encoding="utf-8").write("\n".join(["Brain Vision Data Exchange Marker File, Version 1.0", "", "[Common Infos]", "Codepage=UTF-8", f"DataFile={stem}_ieeg.eeg", "", "[Marker Infos]", "Mk1=New Segment,,1,1,0"]) + "\n")
        typ = []
        for j in range(nch):
            t = etype[j].lower()
            nm = names[j].upper()
            if "EKG" in nm or "ECG" in nm:
                typ.append("ECG")
            elif t == "depth":
                typ.append("SEEG")
            elif t in ("grid", "strip"):
                typ.append("ECOG")
            elif nm in ("", "NAN"):
                typ.append("MISC")
            elif "depth" in longn[j].lower():
                typ.append("SEEG")
            elif "strip" in longn[j].lower() or "grid" in longn[j].lower():
                typ.append("ECOG")
            else:
                typ.append("MISC")  # intracranial contact whose electrode type is not stated in the release
        wtsv(os.path.join(d, stem + "_channels.tsv"),
             ["name", "type", "units", "low_cutoff", "high_cutoff", "sampling_frequency", "status", "status_description", "source_name", "long_name", "electrode_type", "region"],
             [[bnames[j], typ[j], unit, "n/a", "n/a", fs, "good", "n/a", names[j], longn[j], etype[j], region[j]] for j in range(nch)])
        desc = ("Raw iEEG recorded immediately before, during and after continuous 100 Hz direct electrical stimulation" if task == "stim"
                else "Raw iEEG during single-pulse electrical stimulation (before/after continuous stimulation, per the release)")
        par = "; ".join(x for x in [f"amplitude {meta['amplitude_mA']} mA" if meta["amplitude_mA"] else "amplitude not stated in file name",
                                    f"site {meta['site']}" if meta["site"] else "site not stated in file name", f"side {meta['side']}" if meta["side"] else ""] if x)
        wjson(os.path.join(d, stem + "_ieeg.json"), {
            "TaskName": task, "TaskDescription": desc,
            "SamplingFrequency": fs, "PowerLineFrequency": 60, "SoftwareFilters": "n/a", "HardwareFilters": "n/a",
            "Manufacturer": {"natus": "Natus", "tdt": "Tucker-Davis Technologies", "neuroomega": "Alpha Omega (NeuroOmega)"}[meta["system"]],
            "iEEGReference": "n/a (not stated in the release)",
            "RecordingType": "continuous", "RecordingDuration": ns / fs,
            "ECOGChannelCount": typ.count("ECOG"), "SEEGChannelCount": typ.count("SEEG"), "ECGChannelCount": typ.count("ECG"), "MiscChannelCount": typ.count("MISC"),
            "ElectricalStimulation": True,
            "ElectricalStimulationParameters": ("continuous 100 Hz stimulation; " if task == "stim" else "single-pulse stimulation; ") + par + f" (from source file name {bn})",
        })
        wtsv(os.path.join(d, stem + "_events.tsv"), ["onset", "duration", "trial_type"], [[0.0, "%.6f" % (ns / fs), f"recording_{task}"]])
        # in-process round-trip: memmap of the written file vs source dataset (3 windows), mne header read
        import mne
        mm = np.memmap(os.path.join(d, stem + "_ieeg.eeg"), dtype="<f4", mode="r").reshape(ns, nch)
        rt_ok = True
        for s0 in (0, max(0, ns // 2 - 1000), max(0, ns - 2000)):
            a = ds[s0:s0 + 2000, :].astype("<f4")
            rt_ok &= bool(np.array_equal(np.asarray(mm[s0:s0 + 2000]), a, equal_nan=True))
        rr = mne.io.read_raw_brainvision(os.path.join(d, stem + "_ieeg.vhdr"), preload=False, verbose="error")
        rt_ok &= (rr.n_times == ns and len(rr.ch_names) == nch and abs(rr.info["sfreq"] - fs) < 1e-6)
        del mm
        report["runs"].append({"roundtrip_ok": rt_ok, "bids": f"{sub}/ieeg/{stem}_ieeg.vhdr", "source": f"{zn}::{i.filename}", "sha256_source_mat": src_sha, "data_var": dkey,
                               "n_samples": ns, "n_channels": nch, "sfreq": fs, "unit_assumed": unit, "names_from": names_src,
                               "max_abs_float32_rounding_error": maxerr, "max_abs_value": absmax, **meta})
        f.close()
        os.remove(p)
        print("run", stem, ns, nch, fs, maxerr, flush=True)
# electrodes per subject
for sub, info in subjects.items():
    d = os.path.join(OUT, sub, "ieeg")
    el = info["el"]
    if el:
        seen = {}
        rows = []
        for r in el:
            nm = r["name"] or "unnamed"
            seen[nm] = seen.get(nm, 0) + 1
            rows.append([nm if seen[nm] == 1 else f"{nm}_{seen[nm]}", *["%.6f" % v for v in r["xyz"]], "n/a", r["type"] or "n/a", r["long"] or "n/a", r["region"] or "n/a"])
        wtsv(os.path.join(d, f"{sub}_space-Other_electrodes.tsv"), ["name", "x", "y", "z", "size", "type", "long_name", "region"], rows)
        wjson(os.path.join(d, f"{sub}_space-Other_electrodes.json"), {"long_name": {"Description": "Long electrode name (anatomy column 2 of the release electrode file)"},
                                                                     "region": {"Description": "Anatomical region label (anatomy column 4)"}})
        wjson(os.path.join(d, f"{sub}_space-Other_coordsystem.json"), {
            "iEEGCoordinateSystem": "Other", "iEEGCoordinateUnits": "mm",
            "iEEGCoordinateSystemDescription": "elecmatrix of the release electrode file (clinical_elecs_all.mat / TDT_elecs_all.mat), patient-specific space from the authors' MRI/CT reconstruction (README: 'MRI and CT reconstructions allowed us to assign the activity recorded on specific electrodes to different brain regions'); exact space not stated."})
    else:
        names = set()
        for fn in os.listdir(d):
            if fn.endswith("_channels.tsv"):
                names |= {r["name"] for r in csv.DictReader(open(os.path.join(d, fn)), delimiter="\t")}
        wtsv(os.path.join(d, f"{sub}_space-Other_electrodes.tsv"), ["name", "x", "y", "z", "size"], [[n, "n/a", "n/a", "n/a", "n/a"] for n in sorted(names)])
        wjson(os.path.join(d, f"{sub}_space-Other_electrodes.json"), {"SpatialReference": "none: no electrode coordinates in the release for this participant"})
        wjson(os.path.join(d, f"{sub}_space-Other_coordsystem.json"), {"iEEGCoordinateSystem": "Other", "iEEGCoordinateUnits": "n/a",
                                                                      "iEEGCoordinateSystemDescription": "No electrode coordinate file in the release for this participant (x, y, z = n/a)."})
wtsv(os.path.join(OUT, "participants.tsv"), ["participant_id", "age", "sex", "handedness"], [[s, "n/a", "n/a", "n/a"] for s in sorted(subjects)])
os.makedirs(os.path.join(OUT, "code"), exist_ok=True)
shutil.copy(__file__, os.path.join(OUT, "code", os.path.basename(__file__)))
json.dump(report, open(os.path.join(OUT, "code", "conversion_report.json"), "w"), indent=1, default=str)
# sourcedata: everything except the NaturalBehavior zip (absolute timestamps)
dst = os.path.join(OUT, "sourcedata", "dryad-q6vd6wm2")
os.makedirs(dst, exist_ok=True)
for fn in sorted(os.listdir(SRC)):
    pth = os.path.join(SRC, fn)
    if os.path.isfile(pth) and fn != "iEEG_NaturalBehavior.zip" and not fn.endswith(".part"):
        t = os.path.join(dst, fn)
        if not (os.path.exists(t) and os.path.getsize(t) == os.path.getsize(pth)):
            shutil.copyfile(pth, t)
print("DONE", len(report["runs"]), "runs", len(subjects), "subjects", "skipped", report["skipped"])
