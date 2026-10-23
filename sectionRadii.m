% sectionRadii.m
% Reads a curvatureResults_*.mat file and computes, for every section and every
% track, the mean curvature and the corresponding radius (R = 1/k), in three forms:
% R_net (1/|mean k|, left/right turns cancel), R_typical (1/mean|k|) and
% R_median (1/median|k|). Sections are labelled by surface (z>0 or z<0, from the
% filtered z). Curvature is scaled by umPerUnit for physical units. Results are saved
% to sectionRadii_<tag>.xlsx and trackRadii_<tag>.xlsx, where <tag> is taken from the results file name.
% Also plots a histogram of the section mean curvature (meanK), colored by surface.

clear; 
clc;

%% Settings
umPerUnit = 1;     % physical length of one position unit (SET THIS, see note 1)
trimPts   = 0;     % points dropped at each end of every section (0 = none)
minPts    = 5;     % ignore sections/tracks with fewer finite k values

[fn, fp] = uigetfile('*.mat', 'Select curvature results file');
if isequal(fn,0), error('No file selected.'); end
R = load(fullfile(fp,fn));
allTrackResults = R.allTrackResults;
SF = load(R.filteredFile).AllBugs_filt;

%% Collect per-section and per-track statistics
secRows = [];  trkRows = [];   % section: track sec surface n meanK meanAbsK medAbsK reduced
for t = 1:numel(allTrackResults)
    if isempty(allTrackResults{t}) || isempty(allTrackResults{t}.sections), continue; end
    MF = SF{t};  MF = MF(~any(isnan(MF(:,2:4)),2), :);
    secs = allTrackResults{t}.sections;
    kTrack = [];

    for s = 1:numel(secs)
        k = secs(s).k(:);
        if trimPts > 0 && numel(k) > 2*trimPts, k = k(1+trimPts : end-trimPts); end
        k = k(isfinite(k)) / umPerUnit;               % curvature in 1/physical length
        if numel(k) < minPts, continue; end

        inSec = MF(:,1) >= secs(s).startFrame & MF(:,1) <= secs(s).endFrame;
        surf  = 1 + (median(MF(inSec,4)) < 0);        % 1 = z>0, 2 = z<0

        secRows = [secRows; t, s, surf, numel(k), mean(k), mean(abs(k)), ...
                   median(abs(k)), secs(s).stepReduced]; %#ok<AGROW>
        kTrack = [kTrack; k]; %#ok<AGROW>
    end

    if numel(kTrack) >= minPts
        trkRows = [trkRows; t, numel(kTrack), mean(kTrack), mean(abs(kTrack)), ...
                   median(abs(kTrack))]; %#ok<AGROW>
    end
end

%% Build tables (radius = 1/curvature)
S = array2table(secRows, 'VariableNames', {'track','section','surface','nPts', ...
    'meanK','meanAbsK','medAbsK','stepReduced'});
S.R_net     = 1 ./ abs(S.meanK);      % radius of the NET turning (sign cancels)
S.R_typical = 1 ./ S.meanAbsK;        % radius of the typical turning magnitude
S.R_median  = 1 ./ S.medAbsK;         % robust version

T = array2table(trkRows, 'VariableNames', {'track','nPts','meanK','meanAbsK','medAbsK'});
T.R_net     = 1 ./ abs(T.meanK);
T.R_typical = 1 ./ T.meanAbsK;
T.R_median  = 1 ./ T.medAbsK;

%% Histogram of meanK across sections, split by surface
binW  = 0.007;                                         % bin width, in curvature units (1/unit)
lim   = max(0.3, prctile(abs(S.meanK), 99));          % display range, at least +-0.3
edges = -lim : binW : lim;
cols  = {[0.85 0.33 0.10], [0.00 0.45 0.74]};         % surface 1 (z>0), surface 2 (z<0)

figure('Color','w'); hold on;
for g = 1:2
    histogram(S.meanK(S.surface == g), edges, 'FaceColor', cols{g}, ...
        'FaceAlpha', 0.55, 'EdgeColor', 'none');
end
xline(0, 'k:', 'HandleVisibility', 'off');
xlabel('mean curvature of section (1/unit)'); ylabel('number of sections');
legend({'z > 0 surface', 'z < 0 surface'}); grid on;
title('Distribution of section mean curvature');

[cnt, e] = histcounts(S.meanK, edges);
[~, iMax] = max(cnt);
fprintf('Most common meanK bin: %.3f to %.3f (%d of %d sections)\n', ...
    e(iMax), e(iMax+1), cnt(iMax), height(S));

%% Report and save

% Additional summary statistics

fprintf('\n--- Overall section statistics ---\n');

% Typical radius: reciprocal of mean absolute curvature
r = S.R_typical;
r = r(isfinite(r) & r > 0 & r < 1e4);

fprintf('Valid sections: %d\n', numel(r));
fprintf('Median typical radius: %.3f um\n', median(r));
fprintf('Mean typical radius: %.3f um\n', mean(r));
fprintf('Radius IQR (25th-75th): %.3f - %.3f um\n', ...
    prctile(r,25), prctile(r,75));
fprintf('Radius 5th-95th percentiles: %.3f - %.3f um\n', ...
    prctile(r,5), prctile(r,95));

% Absolute curvature: complementary description of turning strength
kAbs = S.meanAbsK;
kAbs = kAbs(isfinite(kAbs) & kAbs > 0);

fprintf('\nMedian mean absolute curvature: %.4f 1/um\n', ...
    median(kAbs));
fprintf('Curvature 25th-75th percentiles: %.4f - %.4f 1/um\n', ...
    prctile(kAbs,25), prctile(kAbs,75));
fprintf('%d sections, %d tracks\n', height(S), height(T));


for g = 1:2
    m = S.surface == g;
    fprintf('Surface %d: %4d sections | median R_typical = %.2f | median R_net = %.2f\n', ...
        g, sum(m), median(S.R_typical(m)), median(S.R_net(m)));
end
[~, base] = fileparts(fn);                              % results file name, no extension
tag = regexprep(base, '^curvatureResults_', '');        % drop the prefix
tag = regexprep(tag, '_\d{8}_\d{6}$', '');              % drop the trailing timestamp
stepTok = regexp(tag, '_step\d+', 'match', 'once');     % e.g. '_step7'
if ~isempty(stepTok)
    tag = strrep(tag, stepTok, '');                     % remove it from the end
    if contains(tag, '_sig')
        tag = regexprep(tag, '_sig', [stepTok '_sig'], 'once');   % put it before sigma
    else
        tag = [tag stepTok];                            % no sigma in name (raw circle fit)
    end
end
writetable(S, ['sectionRadii_' tag '.xlsx']);
writetable(T, ['trackRadii_'   tag '.xlsx']);
