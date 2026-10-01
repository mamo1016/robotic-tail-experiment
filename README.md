# Robotic tail experiment: analysis code

Analysis code for a PhD thesis on the brain and postural responses to covert failures of a
wearable robotic tail (EEG, EMG and foot pressure).

No data are included. Recordings are not publicly available because of the ethics approval.

## Layout

- `codes_improved/main.m` - loads each session, aligns the recordings, filters the EEG,
  cuts trials and computes the ERSP.
- `codes_improved/functions/for_obtainERD/configure_parameters.m` - all settings
  (sampling rates, filters, trial window, baseline). Set `out_dir`, `baseDir` and
  `epochDir` to your own folders.
- `codes_improved/functions/` - processing functions (EEG, ERSP, foot pressure).
- `codes_improved/stats_beta_erd_ttest_N25.m` - participant-level test of the Surprise ERD.
- `codes_improved/scratch/` - trial tables used by the thesis analyses.
- `thesis/analysis/` - figures and analyses added during the thesis corrections
  (per-participant values, processing walkthrough, trial-type prediction, tail video).

Paths are relative: results are written to `output/`, raw recordings are read from `data/`.

Requires MATLAB (Signal Processing and Statistics toolboxes) and Python 3 for the `.py` scripts.
