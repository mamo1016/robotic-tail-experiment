function config=configure_parameters()

% Define all Numbers here.
config = struct();
config.freq_range = 1:1:15;       % Indices for Frequency loop
config.time_range = 1:1:44;       % Indices for Time loop
config.freq_window = 2;           % Width of freq averaging window
config.time_window = 5;           % Width of time averaging window
config.p_threshold = 0.05;        % Significance level
config.r_threshold = 0.45;        % "Strong correlation" threshold
config.remove_outliers = true;    % Toggle outlier removal
config.plot_details = false;      % Toggle detailed debug plots per window

% config.out_dir = "output";  % alternative location
config.out_dir = "output"; % results folder
config.baseDir = 'data'; % Change to your base directory
% config.epochDir = 'output\epochs'; % alternative location
config.epochDir = 'output\epochs'; % % results folder
% Select a subject and session for visualization
config.vis_subject = 1; % Subject to visualize
config.vis_session = 1; % Session to visualize


config.log_file = fopen('eeg_processing_log.txt', 'w');

% Parameters for processing
config.eeg_fs = 1200; % EEG sampling frequency (Hz)
config.emg_fs = 250;  % EMG sampling frequency (Hz)
config.foot_fs = 83;  % Foot sensor sampling frequency (Hz)


% Filter parameters
config.lowcut = 1;    % High-pass cutoff (Hz)
config.highcut = 40;  % Low-pass cutoff (Hz)
config.notch_freq = 50; % Notch filter frequency (Hz)
config.notch_bandwidth = 5; % Notch filter bandwidth (Hz)

% EEG channels of interest (sensorimotor cortex)
config.coi = [13 9 5 15 12 1 14 11 10]; % Channels of interest from file


% Epoch parameters
config.epoch_start = -2;  % Seconds before trigger
config.epoch_end = 4;     % Seconds after trigger
config.baseline_start = -0.5; % Baseline start (seconds)
config.baseline_end = -0.1;   % Baseline end (seconds)

% Beta band parameters
config.beta_range = [13 30];     % Full beta band (Hz)
config.beta_range = [13 19];     % Lower beta band (Hz) - focus of the study
config.analysis_window = [0.5 1.0]; % Analysis window (s) post-movement onset
config.lower_beta = [13 19];     % Lower beta band (Hz) - focus of the study

% Artifact rejection threshold
config.artifact_threshold = 100; % µV

% Check for GPU availability
config.useGPU = (gpuDeviceCount > 0);

config.restDur = 250*20; % do not analyse first x seconds


config.foot_corrdination = [-6; -5; -3; -2; -1; 1; 2; 3; 5; 6];
%% 3. Time-Frequency Decomposition (STFT)
% We calculate the power spectrum for each trial individually
% Parameters for the windowing (tweak these to trade off time vs freq resolution)
config.window_size = round(config.eeg_fs * 0.4);  % 400ms window (standard for low freqs)
config.overlap     = round(config.window_size * 0.95); % 95% overlap for smooth plot
config.nfft        = 2^nextpow2(config.window_size);  % FFT points (power of 2 is faster)
config.baseline_window = [-1.5 -0.5];


% How many times to repeat the resampling (e.g., 500 or 1000)
config.n_boots = 100;


config.foot_fpass = [0.1 60];