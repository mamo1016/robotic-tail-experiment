%% Export LOOCV Consistency Results for AI Analysis
% Purpose: Convert .mat files to human-readable formats (CSV, JSON, TXT)
% so that Claude Opus can inspect and analyse the results

clear; close all; clc;

loocv_dir = 'output/ga_feature_30runs/toResult';
export_dir = 'output/ga_feature_30runs/export_for_ai';

% Create export directory if it doesn't exist
if ~exist(export_dir, 'dir')
    mkdir(export_dir);
end

%% Find latest LOOCV result file
loocv_files = dir(fullfile(loocv_dir, 'results_LOOCV_*.mat'));
if isempty(loocv_files)
    error('No LOOCV result files found in %s', loocv_dir);
end

% Sort by date (newest last)
[~, idx] = sort([loocv_files.datenum]);
loocv_files = loocv_files(idx);
latest_file = fullfile(loocv_dir, loocv_files(end).name);

fprintf('Latest LOOCV file: %s\n', loocv_files(end).name);
fprintf('File timestamp: %s\n', datestr(loocv_files(end).datenum));

%% Load the latest results
load(latest_file, 'results');

%% Part 1: Export structure overview to TXT
txt_file = fullfile(export_dir, '01_STRUCTURE_OVERVIEW.txt');
fid = fopen(txt_file, 'w');

fprintf(fid, '===== LOOCV RESULTS STRUCTURE OVERVIEW =====\n\n');
fprintf(fid, 'File: %s\n', loocv_files(end).name);
fprintf(fid, 'Timestamp: %s\n\n', datestr(loocv_files(end).datenum));

method_names = fieldnames(results);
fprintf(fid, 'Available Methods: %s\n', strjoin(method_names, ', '));
fprintf(fid, 'Number of methods: %d\n\n', length(method_names));

% Check structure of first method
first_method = method_names{1};
first_results = results.(first_method);
fprintf(fid, 'Example structure (method: %s):\n', first_method);
fprintf(fid, '  Type: %s\n', class(first_results));

if isstruct(first_results)
    field_names = fieldnames(first_results);
    fprintf(fid, '  Fields: %s\n', strjoin(field_names, ', '));

    % Get size info
    for f_idx = 1:min(3, length(field_names))
        field = field_names{f_idx};
        data = first_results.(field);
        if isvector(data) || ismatrix(data)
            fprintf(fid, '    - %s: [%s]\n', field, mat2str(size(data)));
        end
    end
end

fprintf(fid, '\n===== TILE INFORMATION =====\n');
fprintf(fid, 'Tiles per phase: 400 (20 freq × 20 time)\n');
fprintf(fid, 'Phases: Control (tiles 1-400), Surprise (401-800), Learning (801-1200)\n');
fprintf(fid, 'Total tiles: 1200\n');

fclose(fid);
fprintf('✓ Exported to: %s\n', txt_file);

%% Part 2: Export method metadata as JSON
json_file = fullfile(export_dir, '02_METADATA.json');

metadata = struct();
metadata.file_name = loocv_files(end).name;
metadata.file_timestamp = datestr(loocv_files(end).datenum);
metadata.methods = method_names;
metadata.num_methods = length(method_names);
metadata.phases = {'Control', 'Surprise', 'Learning'};
metadata.tiles_per_phase = 400;
metadata.total_tiles = 1200;
metadata.freq_bins = 20;
metadata.time_bins = 20;

json_text = jsonencode(metadata);
writelines(json_text, json_file);
fprintf('✓ Exported to: %s\n', json_file);

%% Part 3: Export first method's mask data as CSV (for inspection)
% This shows which tiles were selected in the LOOCV results
csv_file = fullfile(export_dir, '03_SAMPLE_METHOD_MASKS.csv');

% Take first method
sample_method = method_names{1};
sample_data = results.(sample_method);

% Initialize CSV content
csv_content = {};
csv_content{1} = sprintf('Method,Phase,Tile_Index,Selected_X,Selected_Y');

row_idx = 2;
phase_names = {'Control', 'Surprise', 'Learning'};
phase_ranges = {1:400, 401:800, 801:1200};

% Check available mask fields
if isfield(sample_data, 'mask_x')
    mask_x = sample_data.mask_x;
    mask_y = sample_data.mask_y;

    for p = 1:length(phase_names)
        phase_name = phase_names{p};
        tile_range = phase_ranges{p};

        % Get tiles selected in this phase (mask_x and mask_y are binary)
        phase_mask_x = mask_x(tile_range);
        phase_mask_y = mask_y(tile_range);

        % Find tiles where either x or y (or both) are selected
        combined_mask = phase_mask_x | phase_mask_y;
        selected_indices = find(combined_mask == 1);

        fprintf('  %s: %d tiles with selected features out of 400\n', phase_name, length(selected_indices));

        % Export first 20 selected tiles per phase for inspection
        for t_idx = 1:min(20, length(selected_indices))
            tile_num = selected_indices(t_idx);
            x_sel = phase_mask_x(tile_num);
            y_sel = phase_mask_y(tile_num);
            csv_content{row_idx} = sprintf('%s,%s,%d,%d,%d', sample_method, phase_name, tile_num, x_sel, y_sel);
            row_idx = row_idx + 1;
        end
    end

    writecell(csv_content, csv_file);
    fprintf('✓ Exported to: %s\n', csv_file);
else
    fprintf('⚠ Warning: No mask_x/mask_y fields in results structure\n');
    fprintf('  Available fields in %s: %s\n', sample_method, strjoin(fieldnames(sample_data), ', '));
end

%% Part 4: Summary statistics as TXT
summary_file = fullfile(export_dir, '04_SUMMARY_STATISTICS.txt');
fid = fopen(summary_file, 'w');

fprintf(fid, '===== SUMMARY STATISTICS =====\n\n');
fprintf(fid, 'Analysis covers %d methods across 3 phases\n', length(method_names));
fprintf(fid, 'Methods: %s\n\n', strjoin(method_names, ', '));

fprintf(fid, 'Phases and tile ranges:\n');
fprintf(fid, '  Control:  tiles 1-400\n');
fprintf(fid, '  Surprise: tiles 401-800\n');
fprintf(fid, '  Learning: tiles 801-1200\n\n');

fprintf(fid, 'Each tile represents 1 time-frequency bin (20×20 grid per phase)\n');

fclose(fid);
fprintf('✓ Exported to: %s\n', summary_file);

%% Display what was exported
fprintf('\n===== EXPORT COMPLETE =====\n');
fprintf('Files ready for Claude Opus analysis at:\n  %s\n\n', export_dir);
fprintf('Files created:\n');
fprintf('  1. 01_STRUCTURE_OVERVIEW.txt — Structure of results (fields, sizes)\n');
fprintf('  2. 02_METADATA.json — Metadata (methods, phases, dimensions)\n');
fprintf('  3. 03_SAMPLE_METHOD_MASKS.csv — Sample of selected tiles\n');
fprintf('  4. 04_SUMMARY_STATISTICS.txt — Tile grid and phase info\n');
fprintf('\nNext step: Share these files with Claude Opus for analysis\n');
