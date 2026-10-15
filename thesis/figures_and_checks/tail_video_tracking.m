% Tail movement measured from the experiment video (camera behind the participant, 29 fps).
% Examiner: "robot testing data, validation results demonstrating position control performance".
% Method: the median frame shows the tail at rest (hanging to the left). For each swing, at the
% frame where the most pixels have become darker (tail at its right-most position), the centroid
% of the darker pixels gives the peak angle, and the centroid of the pixels that became lighter on
% the left gives the rest angle. Angles are seen from the attachment point: 0 = straight down,
% positive = to the participant's right. Only rows below the hands are used, so arm movement does
% not count. The swing is in the frontal plane, which faces the camera.
clear; clc;
here = fileparts(mfilename('fullpath'));
vid  = fullfile(here, '..', 'lab_media', 'ex_pic', 'ex time stamp.mp4');
fig_dir = fullfile(here, '..', 'revised', 'images', 'results');

v = VideoReader(vid); fps = v.FrameRate;
rows = 470:600; cols = 160:360;            % below the hands
pivot = [250 347];                         % tail attachment (x, y) in frame pixels
G = []; t = [];
while hasFrame(v)
    f = double(readFrame(v))/255; f = 0.299*f(:,:,1) + 0.587*f(:,:,2) + 0.114*f(:,:,3);
    G(:,:,end+1) = f(rows, cols); t(end+1) = v.CurrentTime; %#ok<SAGROW>
end
G = G(:,:,2:end); n = size(G,3); t = t(:);
ref = median(G, 3);
[yy, xx] = ndgrid(rows, cols);
clean = @(m) m & (conv2(double(m), ones(5), 'same') >= 13);
cnt = zeros(n,1);
for k = 1:n, cnt(k) = nnz(clean((ref - G(:,:,k)) > 0.20)); end

% one swing per trial: peaks of the darker-pixel count at least 5 s apart
[~, pk] = findpeaks(cnt, 'MinPeakHeight', 0.5*max(cnt), 'MinPeakDistance', round(5*fps));
S = table();
for i = 1:numel(pk)
    k = pk(i);
    d = clean((ref - G(:,:,k)) > 0.20);
    l = clean((G(:,:,k) - ref) > 0.20) & xx < pivot(1);
    a_pk   = atan2d(mean(xx(d)) - pivot(1), mean(yy(d)) - pivot(2));
    a_rest = atan2d(mean(xx(l)) - pivot(1), mean(yy(l)) - pivot(2));
    S(end+1,:) = table(t(k), a_rest, a_pk, a_pk - a_rest, ...
        'VariableNames', {'peak_time_s','rest_angle_deg','peak_angle_deg','sweep_deg'}); %#ok<SAGROW>
end
writetable(S, fullfile(here, 'tail_video_swings.csv'));
disp(S);
fprintf('Rest %.1f +/- %.1f deg; peak %.1f +/- %.1f deg; sweep %.1f +/- %.1f deg (n = %d swings)\n', ...
    mean(S.rest_angle_deg), std(S.rest_angle_deg), mean(S.peak_angle_deg), std(S.peak_angle_deg), ...
    mean(S.sweep_deg), std(S.sweep_deg), height(S));
save(fullfile(here,'tail_video_tracking.mat'), 't', 'cnt', 'S');
