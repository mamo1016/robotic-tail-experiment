function emg_ave_plot(all_epoch, config)
subs = fieldnames(all_epoch);

emg_combined = struct();
foot_combined = struct();
token_chk = 1;
temp_emg_ex = [];
temp_emg_un = [];
temp_foot_ex = [];
temp_foot_un = [];

for sub_num = 1:length(subs)
    current_name = subs{sub_num};
    tokens = regexp(current_name, 'subject(\d+)', 'tokens');
    sub_id = str2double(tokens{1}{1});
    if sub_id == token_chk
        temp_emg_ex = [temp_emg_ex; all_epoch.(current_name).expected.emg];
        temp_emg_un = [temp_emg_un; all_epoch.(current_name).unexpected.emg];
        temp_foot_ex = [temp_foot_ex; all_epoch.(current_name).expected.foot];
        temp_foot_un = [temp_foot_un; all_epoch.(current_name).unexpected.foot];
    else
        emg_combined.(sprintf('sub_%d', token_chk)).expected = temp_emg_ex;
        emg_combined.(sprintf('sub_%d', token_chk)).unexpected = temp_emg_un;
        foot_combined.(sprintf('sub_%d', token_chk)).expected = temp_foot_ex;
        foot_combined.(sprintf('sub_%d', token_chk)).unexpected = temp_foot_un;
        temp_emg_ex = all_epoch.(current_name).expected.emg;
        temp_emg_un = all_epoch.(current_name).unexpected.emg;
        temp_foot_ex = all_epoch.(current_name).expected.foot;
        temp_foot_un = all_epoch.(current_name).unexpected.foot;
        token_chk = token_chk + 1;
    end
end
emg_combined.(sprintf('sub_%d', token_chk)).expected = temp_emg_ex;
emg_combined.(sprintf('sub_%d', token_chk)).unexpected = temp_emg_un;
foot_combined.(sprintf('sub_%d', token_chk)).expected = temp_foot_ex;
foot_combined.(sprintf('sub_%d', token_chk)).unexpected = temp_foot_un;

sub_names = fieldnames(emg_combined);
min_len = inf;
for s = 1:length(sub_names)
    n = size(emg_combined.(sub_names{s}).unexpected, 1);
    if n < min_len
        min_len = n;
    end
end

n_subs = length(sub_names);
n_samples = size(emg_combined.(sub_names{1}).expected, 2);
str = ["expected", "unexpected"];
emg_ave = zeros(2, n_subs, n_samples);
foot_ave = zeros(2, n_subs, n_samples);

n_boot = 200;
fprintf('Computing bootstrapped EMG/CoP averages (%d iterations)...\n', n_boot);

for s = 1:n_subs
    for c = 1:2
        % EMG bootstrap
        data_emg = emg_combined.(sub_names{s}).(str(c));
        boot_emg = zeros(1, n_samples);
        for b = 1:n_boot
            idx = randi(size(data_emg, 1), [min_len 1]);
            boot_emg = boot_emg + mean(data_emg(idx, :), 1);
        end
        m_emg = boot_emg / n_boot;
        m_emg = m_emg - mean(m_emg(1:1000));
        emg_ave(c, s, :) = m_emg;

        % Foot CoP bootstrap
        data_foot = foot_combined.(sub_names{s}).(str(c));
        boot_foot = zeros(1, n_samples);
        for b = 1:n_boot
            idx_f = randi(size(data_foot, 1), [min_len 1]);
            boot_foot = boot_foot + mean(data_foot(idx_f, :), 1);
        end
        m_foot = boot_foot / n_boot;
        m_foot = m_foot - mean(m_foot(1:1000));
        foot_ave(c, s, :) = m_foot;
    end
end

timeline = config.epoch_start:1/config.eeg_fs:config.epoch_end;
clr = {[0 176/255 240/255], [1 0 0]};
labels = {'Expected Trial', 'Unexpected Trial'};
fig1 = figure('Color', 'w', 'Position', [100 100 1000 600]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
h_emg = gobjects(1, 2);
for c = 1:2
    data = squeeze(emg_ave(c, :, :));
    grand_mean = mean(data, 1) * 5 / 1023;
    sem = std(data, 0, 1) / sqrt(n_subs) * 5 / 1023;
    fill([timeline fliplr(timeline)], [grand_mean + sem, fliplr(grand_mean - sem)], clr{c}, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    h_emg(c) = plot(timeline, grand_mean, 'Color', clr{c}, 'LineWidth', 2);
end
xline(0, 'LineStyle', ':', 'LineWidth', 1, 'HandleVisibility', 'off');
ylabel('EMG Voltage [V]', 'FontSize', 14);
title('Grand-Average EMG', 'FontSize', 16, 'FontWeight', 'bold');
legend(h_emg, labels, 'Location', 'northeast', 'FontSize', 12);
set(gca, 'FontSize', 12);
hold off;

nexttile;
hold on;
h_foot = gobjects(1, 2);
for c = 1:2
    data = squeeze(foot_ave(c, :, :));
    grand_mean = mean(data, 1);
    sem = std(data, 0, 1) / sqrt(n_subs);
    fill([timeline fliplr(timeline)], [grand_mean + sem, fliplr(grand_mean - sem)], clr{c}, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    h_foot(c) = plot(timeline, grand_mean, 'Color', clr{c}, 'LineWidth', 2);
end
xline(0, 'LineStyle', ':', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Time [s]', 'FontSize', 14);
ylabel('CoP Displacement', 'FontSize', 14);
title('Grand-Average Centre of Pressure', 'FontSize', 16, 'FontWeight', 'bold');
legend(h_foot, labels, 'Location', 'northeast', 'FontSize', 12);
set(gca, 'FontSize', 12);
hold off;

file_save = fullfile(config.out_dir, 'all_in_one_pic', 'emg_cop_grand_average.png');
exportgraphics(fig1, file_save, 'Resolution', 300);
fprintf('Saved: %s\n', file_save);
end
