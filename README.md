# Robotic tail experiment: analysis code

Analysis code for a PhD thesis on the brain and postural responses to covert failures of a
wearable robotic tail (EEG, EMG and foot pressure).

No data are included. Recordings are not publicly available because of the ethics approval.

Requires MATLAB (Signal Processing and Statistics toolboxes) and Python 3 for the `.py` scripts.
Paths are relative: run the scripts from the `analysis/` folder. Raw recordings are read from
`data/` and results are written to `output/`. Set `out_dir`, `baseDir` and `epochDir` in
`analysis/functions/preprocessing/configure_parameters.m` to your own folders.

## Order of use

1. `analysis/preprocess_eeg_and_foot.m`: from raw recordings to trials and ERSP (run first).
2. `analysis/test_surprise_erd.m`: the Surprise ERD test reported in the thesis.
3. `analysis/plot_ersp_maps.m`, `analysis/posture_across_session.m`: the main figures.
4. `analysis/trial_tables/`, then `thesis/figures_and_checks/`: figures and analyses added
   during the thesis corrections.

## Folders

| Folder | What it contains |
|---|---|
| `analysis/` | Main scripts (listed below). |
| `analysis/functions/preprocessing/` | Loading the recordings, aligning EEG, EMG and foot pressure, filtering, finding the trigger, cutting trials. Settings in `configure_parameters.m`. |
| `analysis/functions/ersp/` | Time-frequency power (ERSP) of each trial, artefact checks, channel repair, averaging. |
| `analysis/functions/foot_pressure/` | Foot-pressure trials: cleaning, side-to-side position, features, comparisons. |
| `analysis/functions/feature_search/` | Search for brain and foot features that go together (genetic algorithm, tested by leaving one participant out). |
| `analysis/trial_tables/` | Builds one table with a row per trial, used by the thesis foot figures and the trial-type prediction. |
| `thesis/figures_and_checks/` | Figures and checks added during the thesis corrections. |

## Main scripts (`analysis/`)

| Script | What it does |
|---|---|
| `preprocess_eeg_and_foot.m` | Loads each session, aligns the recordings, filters the EEG, cuts trials around the EMG trigger and computes the ERSP and foot-pressure measures. |
| `test_surprise_erd.m` | Surprise ERD (failure minus normal trials, 13–30 Hz, 0–1 s) for each channel, tested across the 25 participants with usable EEG. |
| `test_surprise_erd_26_participants.m` | Earlier version of the same test, including the participant later excluded for poor EEG. |
| `plot_ersp_maps.m` | Time-frequency maps of the brain response. |
| `plot_topomap_surprise_erd.m` | Scalp map of the Surprise ERD. |
| `posture_across_session.m` | How the postural response changes from early to late trials (each participant's trials split into ten equal parts). |
| `posture_across_session_unexpected_trials.m` | Earlier, simpler version of the above for failure trials only. |
| `bayesian_equivalence_learning.m` | Checks whether the change of the postural response across the session is small enough to call "no change". |
| `run_loocv_30times.m` | Runs the brain–foot feature search 30 times. |
| `analyze_loocv_consistency.m` | Which brain features the feature search picks again and again across runs. |
| `analyze_ga_foot_feature_selection.m` | Which foot features the feature search picks, for each method. |
| `analyze_ga_foot_svr_focus.m` | The same, in detail for one method (SVR). |
| `compute_average_loocv_performance.m` | Average result of the 30 feature-search runs. |
| `permutation_test_loocv.m` | Repeats the feature search with shuffled data, to see what results chance alone gives. |
| `permutation_test_loocv_check.m` | Shows how far a running shuffle test has got. |
| `plot_loocv_4methods_trimmed.m`, `plot_loocv_7methods_representative.m`, `plot_loocv_scatter_unified_axes.m`, `plot_loocv_brain_body_single_run.m`, `plot_tile_selection_counts.m` | Figures of the feature-search results. |
| `verify_svr_consensus_tiles.m` | Re-checks the features picked by one method (SVR). |
| `export_consistency_for_analysis.m`, `export_ersp_power_direction.m`, `export_mat_summary.m` | Write feature-search results to text and CSV files for inspection. |

The feature-search scripts belong to a planned paper. They are not used for the results reported
in the thesis.

## Trial tables (`analysis/trial_tables/`)

| Script | What it does |
|---|---|
| `trial_level/build_trial_table.m` | One row per trial: brain response and foot measure, computed as in the main pipeline. |
| `cop_corrected/build_cop_corrected.m` | Recomputes the foot measure as a true centre of pressure (pressure-weighted position). |

## Thesis figures and checks (`thesis/figures_and_checks/`)

| Script | What it does |
|---|---|
| `fig_erd_per_participant.m` | Surprise ERD for each participant; stops if the group values do not match the published ones. |
| `plot_erd_per_participant.m` | Draws that figure from the saved values. |
| `fig_walkthrough.m` | Each processing stage for one typical participant, from EMG to the brain and foot measures. |
| `trial_type_classifier.m` | Can foot pressure tell a failure trial from a normal one? Tested on participants the model has not seen. |
| `tail_video_tracking.m` | Tail angle measured from the experiment video. |
| `tail_video_sheets.m` | Picture sheets of each tail swing, for checking the video measurement by eye. |

Names of the helper functions in `analysis/functions/` describe their task (for example
`extract_foot_epochs.m` cuts foot-pressure trials). They were kept as they are so that the
scripts above still find them.
