% curvatureVsSpeed.m
% Reads a curvatureResults_*.mat file and makes dot plots of mean curvature (x)
% against mean speed (y), one dot per section and one dot per track. Speed is
% computed from the filtered x,y positions (low noise) over each section's frames.
% Sections are coloured by surface (z>0 or z<0, from the filtered z). Curvature
% and speed are scaled by umPerUnit and framesPerSec for physical units.

clear; clc;

%% Settings
umPerUnit    = 1;      % physical length of one position unit (SET THIS)
framesPerSec = 1;      % frame rate (SET THIS), 1 = speed in units per frame
useAbsK      = false;   % true: mean |k| (turning magnitude), false: signed mean k
use3D        = false;  % true: 3D speed (x,y,z), false: in-plane speed (x,y)
minPts       = 5;      % ignore sections with fewer finite k values

[fn, fp] = uigetfile('*.mat', 'Select curvature results file');
if isequal(fn,0), error('No file selected.'); end
R = load(fullfile(fp,fn));
allTrackResults = R.allTrackResults;
SF = load(R.filteredFile).AllBugs_filt;

%% Collect per-section and per-track values
secRows = [];  trkRows = [];    % section: track surface meanK meanV | track: meanK meanV
for t = 1:numel(allTrackResults)
    if isempty(allTrackResults{t}) || isempty(allTrackResults{t}.sections), continue; end
    MF = SF{t};  MF = MF(~any(isnan(MF(:,2:4)),2), :);
    secs = allTrackResults{t}.sections;
    kTrack = [];  vTrack = [];

    for s = 1:numel(secs)
        k = secs(s).k(:);
        k = k(isfinite(k)) / umPerUnit;
        if numel(k) < minPts, continue; end

        inSec = MF(:,1) >= secs(s).startFrame & MF(:,1) <= secs(s).endFrame;
        if sum(inSec) < 2, continue; end
        fr = MF(inSec,1);
        if use3D, P = MF(inSec,2:4); else, P = MF(inSec,2:3); end
        dfr  = diff(fr);
        stp  = sqrt(sum(diff(P).^2, 2));
        ok   = dfr == 1;                                   % never use steps across gaps
        v    = stp(ok) ./ dfr(ok) * umPerUnit * framesPerSec;
        if isempty(v), continue; end

        if useAbsK, kMean = mean(abs(k)); else, kMean = mean(k); end
        surf = 1 + (median(MF(inSec,4)) < 0);              % 1 = z>0, 2 = z<0

        secRows = [secRows; t, surf, kMean, mean(v)]; %#ok<AGROW>
        kTrack = [kTrack; k]; vTrack = [vTrack; v];       %#ok<AGROW>
    end

    if ~isempty(kTrack)
        if useAbsK, kM = mean(abs(kTrack)); else, kM = mean(kTrack); end
        trkRows = [trkRows; t, kM, mean(vTrack)]; %#ok<AGROW>
    end
end

%% Plot: one dot per section, coloured by surface
cols = {[0.85 0.33 0.10], [0.00 0.45 0.74]};
figure('Color','w'); hold on;
for g = 1:2
    m = secRows(:,2) == g;
    scatter(secRows(m,3), secRows(m,4), 18, cols{g}, 'filled', 'MarkerFaceAlpha', 0.5);
end
xlabel('mean curvature of section'); ylabel('mean speed of section');
legend({'z > 0 surface','z < 0 surface'}); grid on;
title('Curvature vs speed (one dot per section)');

%% Plot: one dot per track
figure('Color','w');
scatter(trkRows(:,2), trkRows(:,3), 18, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('mean curvature of track'); ylabel('mean speed of track'); grid on;
title('Curvature vs speed (one dot per track)');

%% Correlation
cS = corrcoef(secRows(:,3), secRows(:,4));
cT = corrcoef(trkRows(:,2), trkRows(:,3));
fprintf('%d sections, Pearson r = %.3f\n', size(secRows,1), cS(1,2));
fprintf('%d tracks,   Pearson r = %.3f\n', size(trkRows,1), cT(1,2));