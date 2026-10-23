% speedHistogram.m
% Reads a curvatureResults_*.mat file and plots two histograms of speed: (1) the mean
% speed of each near-surface section, (2) the single-frame step speeds of all sections
% pooled. Speed is path length per time from the filtered positions over each section's
% frames, in-plane (x,y) by default, same as curvatureVsSpeed.m. Uses the same section
% filter (minPts finite curvature values), so it describes the same sections.
% Prints the mean, median and low percentiles.

clear; clc;

%% Settings
umPerUnit = 1;      % micrometers per position unit
fps       = 15;     % frames per second
use3D     = false;  % true: 3D speed, false: in-plane (x,y)
minPts    = 5;      % same filter as curvatureVsSpeed.m
binW      = 1;      % histogram bin width (um/s)

[fn, fp] = uigetfile('*.mat', 'Select curvature results file');
if isequal(fn,0), error('No file selected.'); end
R = load(fullfile(fp,fn));
allTrackResults = R.allTrackResults;
SF = load(R.filteredFile).AllBugs_filt;

%% Collect speeds
secSpeed = [];  stepSpeed = [];
for t = 1:numel(allTrackResults)
    if isempty(allTrackResults{t}) || isempty(allTrackResults{t}.sections), continue; end
    MF = SF{t};  MF = MF(~any(isnan(MF(:,2:4)),2), :);
    secs = allTrackResults{t}.sections;
    for s = 1:numel(secs)
        k = secs(s).k(:);  k = k(isfinite(k));
        if numel(k) < minPts, continue; end
        inSec = MF(:,1) >= secs(s).startFrame & MF(:,1) <= secs(s).endFrame;
        if sum(inSec) < 2, continue; end
        fr = MF(inSec,1);
        if use3D, P = MF(inSec,2:4); else, P = MF(inSec,2:3); end
        dfr = diff(fr);
        stp = sqrt(sum(diff(P).^2, 2));
        ok  = dfr == 1;
        v   = stp(ok) ./ dfr(ok) * umPerUnit * fps;
        if isempty(v), continue; end
        secSpeed  = [secSpeed;  mean(v)]; %#ok<AGROW>
        stepSpeed = [stepSpeed; v];       %#ok<AGROW>
    end
end

%% Report

% Section-level speeds: each section contributes one value
fprintf('\n--- Section mean speeds ---\n');
fprintf('Number of sections: %d\n', numel(secSpeed));
fprintf('Mean: %.2f um/s\n', mean(secSpeed));
fprintf('Median: %.2f um/s\n', median(secSpeed));
fprintf('Standard deviation: %.2f um/s\n', std(secSpeed));
fprintf('IQR (25th-75th): %.2f - %.2f um/s\n', ...
    prctile(secSpeed,25), prctile(secSpeed,75));
fprintf('Percentiles (5/10/90/95): %.2f %.2f %.2f %.2f um/s\n', ...
    prctile(secSpeed,[5 10 90 95]));

% Modal bin using the same bin width as the histogram
edgesSec = 0:binW:(ceil(max(secSpeed)/binW)*binW + binW);
[countSec, edgesSec] = histcounts(secSpeed, edgesSec);
[~, idxSec] = max(countSec);
fprintf('Modal speed bin: %.2f - %.2f um/s (%d sections)\n', ...
    edgesSec(idxSec), edgesSec(idxSec+1), countSec(idxSec));

% Pooled single-frame speeds: each step contributes one value
fprintf('\n--- Pooled single-frame step speeds ---\n');
fprintf('Number of steps: %d\n', numel(stepSpeed));
fprintf('Mean: %.2f um/s\n', mean(stepSpeed));
fprintf('Median: %.2f um/s\n', median(stepSpeed));
fprintf('Standard deviation: %.2f um/s\n', std(stepSpeed));
fprintf('IQR (25th-75th): %.2f - %.2f um/s\n', ...
    prctile(stepSpeed,25), prctile(stepSpeed,75));
fprintf('Percentiles (5/10/90/95): %.2f %.2f %.2f %.2f um/s\n', ...
    prctile(stepSpeed,[5 10 90 95]));

% Modal bin for pooled step speeds, using the histogram's bin width
edgesStep = 0:binW:(ceil(prctile(stepSpeed,99.5)/binW)*binW + binW);
[countStep, edgesStep] = histcounts(stepSpeed, edgesStep);
[~, idxStep] = max(countStep);
fprintf('Modal speed bin: %.2f - %.2f um/s (%d steps)\n', ...
    edgesStep(idxStep), edgesStep(idxStep+1), countStep(idxStep));

%% Plots
figure('Color','w');
histogram(secSpeed, 0:binW:max(secSpeed)+binW, 'FaceColor', [0.2 0.5 0.8], 'EdgeColor', 'none');
xline(median(secSpeed), 'r-', 'LineWidth', 1.5);
xlabel('mean speed of section (\mum/s)'); ylabel('number of sections'); grid on;
title('Section mean speed (red = median)');

figure('Color','w');
histogram(stepSpeed, 0:binW:prctile(stepSpeed,99.5), 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
xline(median(stepSpeed), 'r-', 'LineWidth', 1.5);
xlabel('single-frame step speed (\mum/s)'); ylabel('number of steps'); grid on;
title('All single-frame step speeds, pooled (red = median)');
