% Contact sheets of the tail region around each swing (0.1 s steps, from -1.6 s to +1.3 s
% relative to the automatically detected peak), for reading the movement start and the
% return to rest by eye. 10 frames per row; row 1 starts at -1.6 s.
here = fileparts(mfilename('fullpath'));
out  = getenv('SHEET_DIR');
T = readtable(fullfile(here,'tail_video_swings.csv'));
v = VideoReader(fullfile(here,'..','lab_media','ex_pic','ex time stamp.mp4'));
for s = 1:height(T)
    tp = T.peak_time_s(s); ts = tp-1.6:0.1:tp+1.3;
    tiles = cell(1, numel(ts));
    for k = 1:numel(ts)
        v.CurrentTime = max(ts(k), 0); f = readFrame(v);
        c = f(340:600, 150:360, :); c(1:4,:,:) = 0; c(:,1:2,:) = 0;   % thin border
        tiles{k} = c;
    end
    nc = 10; nr = ceil(numel(tiles)/nc); [h, w, ~] = size(tiles{1});
    M = uint8(255*ones(nr*h, nc*w, 3));
    for k = 1:numel(tiles)
        r = floor((k-1)/nc); cc = mod(k-1, nc);
        M(r*h+(1:h), cc*w+(1:w), :) = tiles{k};
    end
    ri = round(linspace(1, size(M,1), round(size(M,1)*0.6)));
    ci = round(linspace(1, size(M,2), round(size(M,2)*0.6)));
    imwrite(M(ri, ci, :), fullfile(out, sprintf('sheet_%d.png', s)));
end
