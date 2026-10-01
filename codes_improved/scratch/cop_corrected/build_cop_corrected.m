function build_cop_corrected()
% T2: rebuild the foot measure as a true centre of pressure.
% Plan: paper/notes/COP_PLAN.md. Run 2026-09-09.
%
% Old (published):  load_average = sum(p_i * x_i)          <- a MOMENT
% New (corrected):  cop_norm     = sum(p_i * x_i)/sum(p_i) <- a CoP
%
% Route: reuse the project's own footdist/load_data functions, resample onto the
% epoch file's time_master, and cut at the stored idxs_*. Event detection, EEG
% alignment and the EMG offset are untouched.
%
% VALIDATION GATE: the old moment recomputed through this route must match the
% stored data_resampled_foot. Nothing downstream is trusted unless it does.

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));            % codes_improved
addpath(fullfile(root,'functions','for_obtainERD'));
addpath(fullfile(root,'functions','for_ersp'));
addpath(fullfile(root,'functions','for_foot'));

config = configure_parameters();
ED  = fullfile(config.out_dir,'epochs');
OUT = fullfile(config.out_dir,'trial_level');

fs       = 1200;
t_epoch  = (-2 : 1/fs : 4)';
cop_wins = {[-0.5 0],[0 0.5],[0.5 1.0],[1.0 1.5],[0 1.5]};
cop_names= {'Baseline','Reflex','EarlyComp','LateComp','TotalPost'};
conds    = {'epochs_expected','epochs_unexpected','epochs_right_before','epochs_right_after'};
idxfld   = {'idxs_expected','idxs_unexpected','idxs_right_before','idxs_right_after'};
tcodes   = [1 2 3 4];
KEEP     = 1:25;

subjects = dir(fullfile(config.baseDir,'Sub*'));
rows = {}; val = [];
fprintf('T2: rebuilding CoP as sum(p*x)/sum(p)\n\n');

for s = KEEP
    subject_folder = fullfile(config.baseDir, subjects(s).name);
    files = all_files(subject_folder);
    for ses = 1:length(files.eeg_files)
        fp = fullfile(ED, sprintf('subject_%d___session%d___epoch_data.mat', s, ses));
        if ~exist(fp,'file'), continue; end

        % --- raw foot, decoded by the project's own function ---
        try
            [matchedfiles, ~] = eeg_emg_timeoffset(files, ses, config);
            ds = load_data_optimized(subject_folder, matchedfiles, config);
            ds = footdist(ds);
        catch ME
            fprintf('  sub %d ses %d: SKIP (%s)\n', s, ses, ME.message);
            continue;
        end

        P = double([ds.foot_raw_dist_L, ds.foot_raw_dist_R]);   % [n x 10] pressures
        x = config.foot_corrdination(:);                        % [10 x 1] lateral positions

        moment = P * x;                 % old, published
        total  = sum(P, 2);             % total load
        copn   = moment ./ total;       % new, a true CoP
        copn(total <= 0) = NaN;         % undefined when unloaded

        % same trim as foot_load_average.m
        moment = moment(2:end-1);
        copn   = copn(2:end-1);
        total_t= total(2:end-1);
        tboth  = (ds.foot_time_dist_L + ds.foot_time_dist_R)/2;
        tboth  = tboth(2:end-1);

        E = load(fp); d = E.epoch_data;
        tm = d.time_master;

        old_rs = interp1(tboth, moment,  tm, 'linear','extrap');
        new_rs = interp1(tboth, copn,    tm, 'linear','extrap');
        tot_rs = interp1(tboth, total_t, tm, 'linear','extrap');

        % QC: the foot device stops recording before the EEG does, and interp1's
        % 'extrap' silently invents finite values past that point. Flag it.
        foot_t0 = min(tboth); foot_t1 = max(tboth);
        nlive = sum( (sum(P==0,1)./size(P,1)) <= 0.5 );

        % ---- VALIDATION GATE ----
        ref = d.data_resampled_foot(:);
        g = isfinite(ref) & isfinite(old_rs(:));
        md = max(abs(old_rs(g) - ref(g)));
        rr = corr(old_rs(g), ref(g));
        rel = md / max(1e-12, max(abs(ref(g))));
        val(end+1,:) = [s ses md rr rel]; %#ok<AGROW>
        fprintf('  sub %2d ses %d | validation maxdiff %.3e (rel %.2e) r=%.6f | mean load %.1f\n', ...
                s, ses, md, rel, rr, mean(tot_rs(isfinite(tot_rs))));

        % ---- features per trial, corrected measure ----
        for c = 1:numel(conds)
            if ~isfield(d, conds{c}), continue; end
            ep = d.(conds{c});
            if isempty(ep) || ~isfield(ep,'foot'), continue; end
            idx = d.(idxfld{c});
            nTr = size(ep.foot,1);
            for k = 1:nTr
                a = idx(k,1); b = idx(k,2);
                f  = new_rs(a:b); f = f(:);
                tl = tot_rs(a:b); tl = tl(:);

                r = struct();
                r.subject = s; r.session = ses; r.cond = tcodes(c); r.trial_in_cond = k;
                r.copn_nan = any(isnan(f));
                r.load_mean = mean(tl,'omitnan');

                % QC flags
                r.n_live_sensors = nlive;
                r.frac_extrap    = mean(tm(a:b) < foot_t0 | tm(a:b) > foot_t1);
                r.frac_unloaded  = mean(tl <= 0);
                r.copn_valid     = double(r.frac_extrap == 0 & r.frac_unloaded == 0 & ~r.copn_nan);

                bi = t_epoch >= -0.5 & t_epoch <= 0;
                base_pos = mean(f(bi),'omitnan');
                r.copn_base_pos = base_pos;

                for w = 1:numel(cop_wins)
                    ti = t_epoch >= cop_wins{w}(1) & t_epoch <= cop_wins{w}(2);
                    seg = f(ti);
                    if any(isnan(seg))
                        rv=NaN; av=NaN; dv=NaN; sdv=NaN; pk=NaN;
                    else
                        rv  = sqrt(mean(seg.^2));
                        av  = trapz(abs(seg));
                        dsp = seg - base_pos;
                        dv  = sqrt(mean(dsp.^2));
                        sdv = std(seg);
                        pk  = max(abs(dsp));
                    end
                    r.(sprintf('copn_%s_rms', cop_names{w})) = rv;
                    r.(sprintf('copn_%s_area',cop_names{w})) = av;
                    r.(sprintf('copn_%s_dev', cop_names{w})) = dv;
                    r.(sprintf('copn_%s_sd',  cop_names{w})) = sdv;
                    r.(sprintf('copn_%s_peak',cop_names{w})) = pk;
                end
                rows{end+1} = r; %#ok<AGROW>
            end
        end
        clear ds E d;
    end
end

T = struct2table([rows{:}]);
writetable(T, fullfile(OUT,'cop_corrected.csv'));
V = array2table(val, 'VariableNames', {'subject','session','maxdiff','r','rel'});
writetable(V, fullfile(OUT,'cop_validation.csv'));

fprintf('\n=== VALIDATION GATE ===\n');
fprintf('sessions: %d\n', height(V));
fprintf('correlation old-recomputed vs stored: min %.6f  median %.6f\n', min(V.r), median(V.r));
fprintf('relative max difference:              max %.3e  median %.3e\n', max(V.rel), median(V.rel));
if min(V.r) > 0.9999 && max(V.rel) < 1e-6
    fprintf('GATE PASSED\n');
else
    fprintf('GATE FAILED or MARGINAL - inspect before trusting anything downstream\n');
end
fprintf('\n=== QC ===\n');
fprintf('trials total                          : %d\n', height(T));
fprintf('trials touching extrapolated foot data: %d (%.2f%%)\n', sum(T.frac_extrap>0), 100*mean(T.frac_extrap>0));
fprintf('trials with any unloaded sample       : %d (%.2f%%)\n', sum(T.frac_unloaded>0), 100*mean(T.frac_unloaded>0));
fprintf('trials with NaN CoP                   : %d (%.2f%%)\n', sum(T.copn_nan>0), 100*mean(T.copn_nan>0));
fprintf('trials flagged copn_valid             : %d (%.2f%%)\n', sum(T.copn_valid>0), 100*mean(T.copn_valid>0));
fprintf('live sensors per session: min %d median %d max %d\n', ...
    min(T.n_live_sensors), round(median(T.n_live_sensors)), max(T.n_live_sensors));
fprintf('CoP range check (must lie within the sensor span -6..+6): min %.2f max %.2f\n', ...
    min(T.copn_base_pos), max(T.copn_base_pos));

fprintf('\nwrote %d trials to %s\n', height(T), fullfile(OUT,'cop_corrected.csv'));
end
