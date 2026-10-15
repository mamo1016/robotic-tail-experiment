function fixed_data = channel_interpolation(bad_ch,data)
% bad_ch = [5 9 4]; % Example: Let's say Ch5 (the center) is bad
% [Trials x Samples x Channels] (50-60 x 7201 x 9)
% data = all_epoch.subject01_session01.expected.eeg; 

%% Setup: Define Your Coordinates (Example for a 3x3 Grid)
% Replace these with your actual 10-20 coordinates if needed!
% Format: [X, Y] for channels 1 to 9
chan_coords = [
    1, 1;  1, 2;  1, 3;  % Ch 1-3
    2, 1;  2, 2;  2, 3;  % Ch 4-6
    3, 1;  3, 2;  3, 3   % Ch 7-9
];

%%  Calculate Distances
% num_trial = size(data, 1);
num_sample = size(data, 2);
num_chans = size(data, 3);
distances = zeros(num_chans, 1);

% We create a "Good Channel Mask" to ensure we NEVER use bad data
all_indices = 1:num_chans;
good_ch_indices = setdiff(all_indices, bad_ch);
fixed_data = squeeze(data);
bad_coords = chan_coords(bad_ch, :);

for k = 1:length(bad_ch)
    target_ch = bad_ch(k);
    target_coords = chan_coords(target_ch, :);

    % only loop through the 'good_ch_indices'
    distances = zeros(length(good_ch_indices), 1);

    for i = 1:length(good_ch_indices)
        good_ch = good_ch_indices(i);
        distances(i) = norm(chan_coords(good_ch, :) - target_coords);
    end
    
    %% Step 2: Calculate Weights (Inverse Distance Squared)
    % Closer channels = Smaller distance = Higher Weight
    weights = 1 ./ (distances .^ 2);
    
    % Normalize weights so they sum to 1
    weights = weights / sum(weights); 
    
    
    %% Step 3: Apply Interpolation to Data
    % We can vectorize this to do all timepoints/trials at once
    % Formula: New_Ch = Sum(Neighbor_Data * Neighbor_Weight)
    
    % Initialize reconstructed channel
    reconstructed_data = zeros(1, num_sample); 

    for i = 1:length(good_ch_indices)
        good_ch = good_ch_indices(i);

        % Add this channel's contribution: Data * Weight
        neighbor_data = squeeze(data(:, :, good_ch)); 
        reconstructed_data = reconstructed_data + (neighbor_data * weights(i));
    end
    % Store the fixed signal into our new matrix
    fixed_data(:, target_ch) = reconstructed_data;
end