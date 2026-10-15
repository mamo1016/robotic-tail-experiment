function filter_para=design_filter(config)
% Design filters once (outside the loop for efficiency)
filter_para=struct();
[filter_para.b_high, filter_para.a_high] = butter(4, config.lowcut/(config.eeg_fs/2), 'high');
[filter_para.b_low, filter_para.a_low] = butter(6, config.highcut/(config.eeg_fs/2), 'low');
[filter_para.b_notch, filter_para.a_notch] = iirnotch(config.notch_freq/(config.eeg_fs/2), config.notch_freq/(config.notch_bandwidth*config.eeg_fs/2));