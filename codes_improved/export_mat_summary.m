mat_file =
    'output/ga_feature/results_LOOCV_2026-02-24_12-49.mat';

if
  ~exist(mat_file, 'file') fprintf('Error: File not found: %s\n', mat_file);
exit;
end

    data = load(mat_file);

fprintf('=== MAT FILE SUMMARY ===\n');
fprintf('Timestamp: %s\n', data.run_info.timestamp);
fprintf('Total Subjects: %d\n',
        length(data.results.(data.run_info.methods{1}).actual));

methods = data.run_info.methods;
for
  i = 1 : length(methods) m = methods{i};
actual = data.results.(m).actual;
predicted = data.results.(m).predicted;

% Calculate correlation[r, p] =
    corr(predicted, actual, 'Type', 'Spearman', 'Rows', 'complete');

% Count selected features mask_x = data.results.(m).mask_x;
mask_y = data.results.(m).mask_y;

fprintf('\nMethod: %s\n', upper(m));
fprintf('  Correlation: r = %.4f (p = %.4f)\n', r, p);
fprintf('  Selected Brain (EEG) Features: %d / %d\n', sum(mask_x),
        length(mask_x));
fprintf('  Selected Foot Features: %d / %d\n', sum(mask_y), length(mask_y));
end fprintf('========================\n');
