function build_trial_table()
% Build a single-trial table from the raw epoch files.
% ERD replicates functions/for_ersp/ersp.m exactly (per-trial baseline, dB).
% CoP replicates functions/for_foot/generate_foot_features.m minus the trial average.

ED   = 'output\epochs';
OUT  = 'output\trial_level';
if ~exist(OUT,'dir'), mkdir(OUT); end

fs        = 1200;
t_epoch   = (-2 : 1/fs : 4)';            % 7201
win_size  = round(fs*0.4);               % 480
noverlap  = round(win_size*0.95);        % 456
nfft      = 2^nextpow2(win_size);        % 512
base_win  = [-1.5 -0.5];
beta_rng  = [13 30];
erd_win   = [0 1.0];
alpha_rng = [8 13];                      % exploratory only
lbeta_rng = [13 19];                     % exploratory only

cop_wins  = {[-0.5 0],[0 0.5],[0.5 1.0],[1.0 1.5],[0 1.5]};
cop_names = {'Baseline','Reflex','EarlyComp','LateComp','TotalPost'};
conds     = {'epochs_expected','epochs_unexpected','epochs_right_before','epochs_right_after'};
tcodes    = [1 2 3 4];

KEEP = 1:25;                             % participant 26 excluded (poor EEG)

rows = {};
t0 = tic;
for s = KEEP
    for ses = 1:3
        fp = fullfile(ED, sprintf('subject_%d___session%d___epoch_data.mat', s, ses));
        if ~exist(fp,'file'), continue; end
        E = load(fp); d = E.epoch_data;

        for c = 1:numel(conds)
            if ~isfield(d, conds{c}), continue; end
            ep = d.(conds{c});
            if isempty(ep) || ~isfield(ep,'eeg'), continue; end
            nTr = size(ep.eeg,1);
            % onset times for this condition, for global ordering
            tf = strrep(conds{c},'epochs_','times_eeg_st_');
            if isfield(d, tf), onset = d.(tf); else, onset = nan(nTr,1); end

            for k = 1:nTr
                r = struct();
                r.subject = s; r.session = ses; r.cond = tcodes(c);
                r.trial_in_cond = k;
                r.onset = onset(min(k,numel(onset)));

                % ---------- EEG ----------
                eeg = squeeze(ep.eeg(k,:,:));      % [7201 x 9]
                seg = t_epoch >= -1.5 & t_epoch <= 1.5;
                r.eeg_absmax = max(max(abs(eeg(seg,:))));
                r.eeg_nan    = any(any(isnan(eeg(seg,:))));

                bt = nan(1,9); at = nan(1,9); lb = nan(1,9);
                for ch = 1:9
                    x = double(eeg(:,ch));
                    if any(isnan(x)), continue; end
                    [S,F,T] = spectrogram(x, win_size, noverlap, nfft, fs);
                    P  = abs(S).^2;
                    tr = T - 2;                    % align to epoch start of -2 s
                    bi = tr >= base_win(1) & tr <= base_win(2);
                    if ~any(bi), continue; end
                    bv = mean(P(:,bi), 2);
                    db = 10*log10(P ./ repmat(bv,1,size(P,2)));
                    wi = tr >= erd_win(1) & tr <= erd_win(2);
                    fb = F >= beta_rng(1)  & F <= beta_rng(2);
                    fa = F >= alpha_rng(1) & F <= alpha_rng(2);
                    fl = F >= lbeta_rng(1) & F <= lbeta_rng(2);
                    bt(ch) = mean(mean(db(fb,wi)));
                    at(ch) = mean(mean(db(fa,wi)));
                    lb(ch) = mean(mean(db(fl,wi)));
                end
                for ch = 1:9
                    r.(sprintf('erd_ch%d',ch))    = bt(ch);
                    r.(sprintf('alpha_ch%d',ch))  = at(ch);
                    r.(sprintf('lbeta_ch%d',ch))  = lb(ch);
                end

                % ---------- CoP ----------
                f = double(ep.foot(k,:))';
                r.foot_nan = any(isnan(f));
                bi_cop = t_epoch >= -0.5 & t_epoch <= 0;
                base_pos = mean(f(bi_cop),'omitnan');   % pre-perturbation standing position
                r.cop_base_pos = base_pos;
                for w = 1:numel(cop_wins)
                    ti = t_epoch >= cop_wins{w}(1) & t_epoch <= cop_wins{w}(2);
                    seg2 = f(ti);
                    if any(isnan(seg2))
                        rv=NaN; av=NaN; dv=NaN; sdv=NaN; pk=NaN;
                    else
                        rv  = sqrt(mean(seg2.^2));            % as published (DC dominated)
                        av  = trapz(abs(seg2));               % as published
                        dsp = seg2 - base_pos;               % displacement from standing position
                        dv  = sqrt(mean(dsp.^2));
                        sdv = std(seg2);
                        pk  = max(abs(dsp));
                    end
                    r.(sprintf('cop_%s_rms',  cop_names{w})) = rv;
                    r.(sprintf('cop_%s_area', cop_names{w})) = av;
                    r.(sprintf('cop_%s_dev',  cop_names{w})) = dv;
                    r.(sprintf('cop_%s_sd',   cop_names{w})) = sdv;
                    r.(sprintf('cop_%s_peak', cop_names{w})) = pk;
                end

                % ---------- EMG / tail ----------
                if isfield(ep,'emg')
                    g = double(ep.emg(k,:))';
                    for w = 1:numel(cop_wins)
                        ti = t_epoch >= cop_wins{w}(1) & t_epoch <= cop_wins{w}(2);
                        r.(sprintf('emg_%s_rms',cop_names{w})) = sqrt(mean(g(ti).^2));
                    end
                end
                if isfield(ep,'tail')
                    tl = double(ep.tail(k,:))';
                    ti = t_epoch >= 0 & t_epoch <= 1.5;
                    r.tail_post_rms = sqrt(mean(tl(ti).^2));
                end

                rows{end+1} = r; %#ok<AGROW>
            end
        end
        fprintf('sub %2d ses %d done (%d rows, %.1f min)\n', s, ses, numel(rows), toc(t0)/60);
    end
end

T = struct2table([rows{:}]);
% global trial order within session, by onset time
T.trial_overall = nan(height(T),1);
key = T.subject*10 + T.session;
for u = unique(key)'
    idx = find(key==u);
    [~,ord] = sort(T.onset(idx));
    rank = nan(numel(idx),1); rank(ord) = 1:numel(idx);
    T.trial_overall(idx) = rank;
end

save(fullfile(OUT,'trial_table.mat'),'T','-v7.3');
writetable(T, fullfile(OUT,'trial_table.csv'));
fprintf('\nSAVED %d rows, %d columns -> %s\n', height(T), width(T), OUT);
fprintf('rows per condition: %s\n', mat2str(histcounts(T.cond,0.5:4.5)));
end
