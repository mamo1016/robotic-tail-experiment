function chk_dir(output_folder)
    % output_folder = "output\interpolate"
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
        fprintf('Folder created: %s\n', output_folder);
    end