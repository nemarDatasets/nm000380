# Stimulation onsets detected from the stimulation artefact (nm000380)

**These events are NOT from a source.** The release does not annotate them. They were detected in Oct 2026 by the NEMAR curators from the stimulation artefact visible in the recorded signals, after an external review (J. Gadonneix, Oct 2026) noted that the artefact places them. They are an analysis product: check them against your own artefact detection before relying on them. The raw dataset's own `events.tsv` files are unchanged.

## Method

Code: `code/fix-20261030_stimartefact/` of the raw dataset (`stimfeat.py` features, `stimdetect.py` rules, `build_derivative.py` this tree). Data were read by annex key from the NEMAR S3 bucket (size and SHA-256 checked).

1. Contacts: channels typed ECOG, SEEG or DBS in `channels.tsv` with status not bad (recordings whose contacts are all typed MISC use the MISC channels).
2. Step score per sample: 90th percentile across contacts of |x[n+1] - x[n]| divided by that contact's robust SD of steps (1.4826 x median absolute deviation). Stimulation pulses are steps far above the normal sample-to-sample change in most contacts at once.
3. Single pulses: clusters of samples with step score > 50, at least 1 s (SPES) or 0.3 s (SPS) apart. A candidate is kept only if it is part of a chain of at least 5 pulses spaced at the protocol interval stated in the article (SPES: 15 s +- 0.5 s; SPS: 1 s +- 0.1 s). Other candidates are rejected and listed below.
4. 43.3 Hz trains (nm000372 task-hfs): step score > 50 at least 10 ms apart, grouped while gaps are <= 0.1 s; a group with >= 20 pulses and a median interval of 23.1 +- 2 ms (43.3 Hz, article) is a train.
5. 100 Hz blocks (nm000380 task-stim): signals resampled to about 1 kHz; spectrogram (1-s Hann windows, 0.1-s steps), median across contacts of log10 power; comb score = mean over 100, 200 and 400 Hz of log10 P(f) minus the mean log10 power 6-10 Hz away. Frames with comb score > 0.5 (power at the stimulation lines more than 3.2 x the neighbouring bins) count as stimulation; gaps <= 5 s are bridged and blocks shorter than 20 s ignored. Onset and offset are then refined to the first/last step (score > 20, >= 4 ms apart) followed/preceded by >= 8 steps within 100 ms, searched within 1.5 s of the spectral edge (`method` = pulse_refined); otherwise the spectral edge is used (`method` = spectral, about +-0.5 s).

Thresholds were set once from the score distributions of these recordings (background step score: median 0.7-1.4) and are the same for every run.

## Results per recording

- `sub-EC105_task-stim_acq-natus_run-01`: 1 block(s) onset 29.521 s, duration 100.7 s (pulse_refined). Comb score median 0.12, 95th percentile 1.93.
- `sub-EC105_task-stim_acq-natus_run-02`: 1 block(s) onset 27.871 s, duration 104.7 s (pulse_refined). Comb score median 0.16, 95th percentile 1.39.
- `sub-EC105_task-stim_acq-natus_run-03`: 1 block(s) onset 23.359 s, duration 105.2 s (pulse_refined). Comb score median 0.20, 95th percentile 2.97.
- `sub-EC105_task-stim_acq-natus_run-04`: 1 block(s) onset 31.105 s, duration 109.9 s (pulse_refined). Comb score median 1.75, 95th percentile 3.62.
- `sub-EC105_task-stim_acq-natus_run-05`: 1 block(s) onset 21.843 s, duration 104.4 s (pulse_refined). Comb score median 2.65, 95th percentile 3.17.
- `sub-EC108_task-stim_acq-natus_run-01`: 4 block(s) onset 57.777 s, duration 92.1 s (spectral); onset 152.306 s, duration 28.2 s (pulse_refined); onset 182.787 s, duration 47.5 s (pulse_refined); onset 237.173 s, duration 42.2 s (pulse_refined). Comb score median 0.19, 95th percentile 0.83.
- `sub-EC108_task-stim_acq-natus_run-02`: 1 block(s) onset 65.671 s, duration 451.9 s (spectral). Comb score median 1.95, 95th percentile 2.28.
- `sub-EC122_task-stim_acq-natus_run-01`: 1 block(s) onset 27.159 s, duration 113.8 s (pulse_refined). Comb score median 0.11, 95th percentile 3.75.
- `sub-EC125_task-stim_acq-natus_run-01`: 1 block(s) onset 138.462 s, duration 239.5 s (pulse_refined). Comb score median 3.70, 95th percentile 4.18.
- `sub-EC125_task-stim_acq-natus_run-02`: 1 block(s) onset 25.625 s, duration 234.5 s (pulse_refined). Comb score median 2.35, 95th percentile 4.52.
- `sub-EC125_task-stim_acq-natus_run-03`: 1 block(s) onset 37.695 s, duration 190.9 s (pulse_refined). Comb score median 2.45, 95th percentile 2.58.
- `sub-EC129_task-stim_acq-natus_run-01`: 1 block(s) onset 29.058 s, duration 276.3 s (pulse_refined). Comb score median 2.58, 95th percentile 2.72.
- `sub-EC129_task-stim_acq-natus_run-02`: 1 block(s) onset 30.474 s, duration 213.6 s (pulse_refined). Comb score median 2.59, 95th percentile 2.76.
- `sub-EC129_task-stim_acq-natus_run-03`: 1 block(s) onset 44.183 s, duration 211.0 s (pulse_refined). Comb score median 2.82, 95th percentile 3.03.
- `sub-EC133_task-stim_acq-natus_run-01`: 1 block(s) onset 24.068 s, duration 188.0 s (pulse_refined). Comb score median 2.47, 95th percentile 2.60.
- `sub-EC133_task-stim_acq-natus_run-02`: 1 block(s) onset 14.813 s, duration 163.2 s (pulse_refined). Comb score median 3.73, 95th percentile 4.09.
- `sub-EC133_task-stim_acq-natus_run-03`: 1 block(s) onset 18.046 s, duration 157.9 s (pulse_refined). Comb score median 2.50, 95th percentile 2.79.
- `sub-EC133_task-stim_acq-natus_run-04`: 1 block(s) onset 13.584 s, duration 146.5 s (pulse_refined). Comb score median 0.12, 95th percentile 2.32.
- `sub-EC137_task-sps_acq-natus_run-01`: 40 pulses kept of 41 candidates, in 2 chain(s) of 20, 20; rejected candidate onsets (s): 97.188.
- `sub-EC137_task-stim_acq-natus_run-01`: 1 block(s) onset 32.312 s, duration 100.8 s (pulse_refined). Comb score median 1.52, 95th percentile 1.73.
- `sub-EC137_task-stim_acq-natus_run-02`: 1 block(s) onset 53.594 s, duration 100.5 s (pulse_refined). Comb score median 4.04, 95th percentile 4.30.
- `sub-EC137_task-stim_acq-natus_run-03`: 1 block(s) onset 33.799 s, duration 100.5 s (pulse_refined). Comb score median 4.12, 95th percentile 4.36.
- `sub-EC137_task-stim_acq-natus_run-04`: 1 block(s) onset 32.831 s, duration 103.1 s (spectral). Comb score median 4.11, 95th percentile 4.33.
- `sub-EC139_task-sps_acq-natus_run-01`: 40 pulses kept of 41 candidates, in 2 chain(s) of 20, 20; rejected candidate onsets (s): 220.265.
- `sub-EC142_task-stim_acq-tdt_run-01`: 1 block(s) onset 21.475 s, duration 140.5 s (pulse_refined). Comb score median 5.37, 95th percentile 5.50.
- `sub-EC142_task-stim_acq-tdt_run-02`: 1 block(s) onset 16.075 s, duration 124.3 s (pulse_refined). Comb score median 5.14, 95th percentile 5.37.
- `sub-EC142_task-stim_acq-tdt_run-03`: 1 block(s) onset 20.042 s, duration 96.6 s (pulse_refined). Comb score median 3.99, 95th percentile 4.10.
- `sub-EC142_task-stim_acq-tdt_run-04`: 1 block(s) onset 12.213 s, duration 101.8 s (pulse_refined). Comb score median 4.66, 95th percentile 4.98.
- `sub-EC142_task-stim_acq-tdt_run-05`: 1 block(s) onset 9.072 s, duration 90.8 s (pulse_refined). Comb score median 5.00, 95th percentile 5.24.
- `sub-EC150_task-sps_acq-natus_run-01`: 40 pulses kept of 41 candidates, in 2 chain(s) of 20, 20; rejected candidate onsets (s): 140.867.
- `sub-EC150_task-stim_acq-neuroomega_run-01`: 1 block(s) onset 36.240 s, duration 199.2 s (pulse_refined). Comb score median 3.67, 95th percentile 3.76.
- `sub-EC150_task-stim_acq-neuroomega_run-02`: 1 block(s) onset 34.908 s, duration 210.9 s (pulse_refined). Comb score median 3.43, 95th percentile 3.56.
- `sub-EC150_task-stim_acq-neuroomega_run-03`: 1 block(s) onset 44.888 s, duration 215.5 s (pulse_refined). Comb score median 0.34, 95th percentile 3.48.
- `sub-EC151_task-stim_acq-natus_run-01`: 1 block(s) onset 28.844 s, duration 225.5 s (spectral). Comb score median 0.75, 95th percentile 1.06.
- `sub-EC152_task-stim_acq-neuroomega_run-01`: 1 block(s) onset 214.332 s, duration 173.6 s (pulse_refined). Comb score median 0.25, 95th percentile 1.31.
- `sub-EC153_task-sps_acq-natus_run-01`: 120 pulses kept of 133 candidates, in 6 chain(s) of 20, 20, 20, 20, 20, 20; rejected candidate onsets (s): 0.819, 1.496, 264.256, 414.681, 741.407, 921.016, 1070.525, 1201.084, 1210.946, 1242.127, 1276.889, 1303.794, 1346.414.
- `sub-EC153_task-stim_acq-natus_run-01`: 1 block(s) onset 52.797 s, duration 30.7 s (pulse_refined). Comb score median 0.07, 95th percentile 3.03.
- `sub-EC153_task-stim_acq-natus_run-02`: 1 block(s) onset 15.932 s, duration 29.2 s (pulse_refined). Comb score median 0.19, 95th percentile 3.77.
- `sub-EC153_task-stim_acq-natus_run-03`: 1 block(s) onset 24.848 s, duration 38.1 s (pulse_refined). Comb score median 0.20, 95th percentile 2.86.
- `sub-EC153_task-stim_acq-natus_run-04`: 1 block(s) onset 64.829 s, duration 33.1 s (pulse_refined). Comb score median 0.08, 95th percentile 3.94.
- `sub-EC153_task-stim_acq-natus_run-05`: 1 block(s) onset 10.904 s, duration 39.9 s (pulse_refined). Comb score median 0.11, 95th percentile 3.78.
- `sub-EC155_task-sps_acq-natus_run-01`: 40 pulses kept of 43 candidates, in 2 chain(s) of 20, 20; rejected candidate onsets (s): 47.177, 88.278, 302.793.
- `sub-EC162_task-stim_acq-natus_run-01`: 1 block(s) onset 28.296 s, duration 182.8 s (pulse_refined). Comb score median 2.20, 95th percentile 2.39.
- `sub-EC175_task-stim_acq-tdt_run-01`: 1 block(s) onset 24.233 s, duration 254.8 s (pulse_refined). Comb score median 1.72, 95th percentile 2.16.
- `sub-EC99_task-stim_acq-natus_run-01`: 1 block(s) onset 43.315 s, duration 222.5 s (pulse_refined). Comb score median 2.48, 95th percentile 2.71.

## Notes

- task-sps: 280 pulses in 14 trains of 20 pulses at 1 s intervals over the 5 single-pulse recordings (sub-EC137, EC139, EC150 and EC155: 2 trains each; sub-EC153: 6 trains), matching the article's protocol (1 Hz, 20 s). The external review counted 276 of 280.
- task-stim: one 100 Hz block in 39 of the 40 continuous-stimulation recordings. In sub-EC108 run-01 (1 mA) the artefact is weak (comb score about 0.2-0.8 from about 55 s to 279 s, threshold 0.5), and the fixed rules split it into 4 segments; it is most likely one block from about 57.8 s to about 279 s, but the table keeps the 4 detected segments as detected.
- The sham (0 mA) block that the article describes at the start of each session leaves no artefact and cannot be detected this way.
- Block durations range from about 29 s (sub-EC153) to about 452 s (sub-EC108 run-02); the article states 100-200 s for most blocks. Onset/offset of each block are in the table; the stimulation site and amplitude of each run are in the raw sidecar (`ElectricalStimulationParameters`).
