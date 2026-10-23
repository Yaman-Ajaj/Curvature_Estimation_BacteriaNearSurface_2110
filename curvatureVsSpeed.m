% curvatureVsSpeed.m
% Reads a curvatureResults_*.mat file and makes a dot plot of mean curvature (x)
% against mean speed (y), one dot per section. Speed is computed from the filtered
% x,y positions (low noise) over each section's frames. Optionally, dots are coloured
% by section length. Curvature and speed are scaled by umPerUnit and framesPerSec.
%
% Tests whether curvature falls with speed (hypothesis k = omega/v, i.e. b = -1):
%   - power-law fit  log k = log A + b log v
%     MAIN estimator: robust fit (robustfit), with a by-track bootstrap 95% CI
%     COMPARISON:     ordinary least squares (polyfit)
%   - Spearman rho, binned medians, slowest-1% comparison, clean-sections-only refit.
% With cleanOnly = true, sections shorter than minFrames or with a reduced
% analysis step are excluded from everything.
% Saves speeds and lengths of the analysed sections to realSections.mat
% (used by synthCurvatureSpeedTest.m). Requires the Statistics Toolbox.

clear;
rng(1);   % reproducible bootstrap

%% Settings
umPerUnit     = 1;      % physical length of one position unit (SET THIS)
framesPerSec  = 15;     % frame rate (SET THIS), 1 = speed in units per frame
useAbsK       = false;  % true: mean |k| (turning magnitude), false: signed mean k
use3D         = false;  % true: 3D speed (x,y,z), false: in-plane speed (x,y)
minPts        = 5;      % ignore sections with fewer finite k values
xLimit        = 2;      % x-axis limit: [-xLimit xLimit], or [0 xLimit] if useAbsK
colorByLength = true;   % true: colour dots by section length + colour bar
cleanOnly     = false;   % true: only sections with nFrames >= minFrames and no reduced step
minFrames     = 15;
nBoot         = 1000;   % bootstrap resamples (by track)

[fn, fp] = uigetfile('*.mat', 'Select curvature results file');
if isequal(fn,0), error('No file selected.'); end
R = load(fullfile(fp,fn));
allTrackResults = R.allTrackResults;
SF = load(R.filteredFile).AllBugs_filt;

%% Collect per-section values: [track, meanK, meanV, nFrames, stepReduced, surface]
secRows = [];
for t = 1:numel(allTrackResults)
    if isempty(allTrackResults{t}) || isempty(allTrackResults{t}.sections), continue; end
    MF = SF{t};  MF = MF(~any(isnan(MF(:,2:4)),2), :);
    secs = allTrackResults{t}.sections;

    for s = 1:numel(secs)
        k = secs(s).k(:);
        k = k(isfinite(k)) / umPerUnit;
        if numel(k) < minPts, continue; end

        inSec = MF(:,1) >= secs(s).startFrame & MF(:,1) <= secs(s).endFrame;
        nIn = sum(inSec);
        if nIn < 2, continue; end
        if cleanOnly && (secs(s).stepReduced || nIn < minFrames), continue; end

        fr = MF(inSec,1);
        if use3D, P = MF(inSec,2:4); else, P = MF(inSec,2:3); end
        dfr = diff(fr);
        stp = sqrt(sum(diff(P).^2, 2));
        ok  = dfr == 1;                                   % never use steps across gaps
        v   = stp(ok) ./ dfr(ok) * umPerUnit * framesPerSec;
        if isempty(v), continue; end

        if useAbsK, kMean = mean(abs(k)); else, kMean = mean(k); end
        surf = 1 + (median(MF(inSec,4)) < 0);             % 1 = z>0, 2 = z<0
        secRows = [secRows; t, kMean, mean(v), nIn, secs(s).stepReduced, surf]; %#ok<AGROW>
    end
end

%% Save real speeds and section lengths for the synthetic test
vReal   = secRows(:,3) / framesPerSec / umPerUnit;   % position units per frame
nFrReal = secRows(:,4);                              % section length in frames
save('realSections.mat', 'vReal', 'nFrReal');

%% Dot plot (x-limit applied to the plot and Pearson r only)
if useAbsK, xRange = [0 xLimit]; else, xRange = [-xLimit xLimit]; end
shown = secRows(:,2) >= xRange(1) & secRows(:,2) <= xRange(2);
fprintf('%d sections, %d outside the x-limit and not shown\n', size(secRows,1), sum(~shown));
D = secRows(shown,:);

figure('Color','w'); hold on;
if colorByLength
    scatter(D(:,2), D(:,3), 18, D(:,4), 'filled', 'MarkerFaceAlpha', 0.6);
    colormap(parula);
    cb = colorbar;  cb.Label.String = 'section length (frames)';
    clim([0 150]);
else
    scatter(D(:,2), D(:,3), 18, [0.2 0.4 0.7], 'filled', 'MarkerFaceAlpha', 0.5);
end
xlim(xRange);
xlabel('mean curvature of section'); ylabel('mean speed of section'); grid on;
title('Curvature vs speed (one dot per section)');

c = corrcoef(D(:,2), D(:,3));
fprintf('Pearson r = %.3f (N = %d)\n', c(1,2), size(D,1));

%% Prepare data for the power-law analysis (all sections, no x-limit)
k = abs(secRows(:,2));  v = secRows(:,3);
ok = k > 0 & v > 0 & isfinite(k) & isfinite(v);
k = k(ok);  v = v(ok);
trk = secRows(ok,1);  nFr = secRows(ok,4);  red = secRows(ok,5) == 1;  sf = secRows(ok,6);

[~, ~, g] = unique(trk);
nPer = accumarray(g, 1);
fprintf('sections per track: median %d, mean %.1f, max %d\n', median(nPer), mean(nPer), max(nPer));
fprintf('speed percentiles 1/2/5/10: %.2f %.2f %.2f %.2f\n', prctile(v,[1 2 5 10]));

%% Power-law fits: robust (main) and OLS (comparison); b = -1 means k = omega/v
fprintf('\n--- All sections ---\n');
pAll = fitReport(v, k, trk, 'all sections', nBoot);
for q = [1 2 5 10]                                    % does b depend on the slowest sections?
    m = v >= prctile(v, q);
    fitReport(v(m), k(m), trk(m), sprintf('drop slowest %d%%', q), 0);
end

%% Median curvature in speed bins (10 bins, equal numbers of sections)
edges = prctile(v, 0:10:100);
[~, ~, bin] = histcounts(v, edges);
vMid = nan(10,1);  kMed = nan(10,1);
for i = 1:10
    vMid(i) = median(v(bin == i));
    kMed(i) = median(k(bin == i));
end
slopeBetweenBins = diff(log(kMed)) ./ diff(log(vMid));
disp([vMid(1:end-1), vMid(2:end), slopeBetweenBins]);   % local exponent in each speed range

%% Log-log plot
vv = linspace(min(v), max(v), 100);
figure('Color','w');
loglog(v, k, '.', 'Color', [0.7 0.7 0.7]); hold on;
loglog(vMid, kMed, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
loglog(vv, exp(pAll(1)) * vv.^pAll(2), 'r-', 'LineWidth', 1.5);
pOLS = polyfit(log(v), log(k), 1);
loglog(vv, exp(pOLS(2)) * vv.^pOLS(1), 'b--', 'LineWidth', 1);
xlabel('mean speed of section'); ylabel('mean |curvature| of section'); grid on;
legend({'sections','median per speed bin','robust fit','OLS fit'}, 'Location', 'southwest');
title('Curvature vs speed (log-log)');

%% Who are the slowest 1%?
slow = v <= prctile(v, 1);
fprintf('\n%-22s %10s %10s\n', '', 'slowest 1%', 'the rest');
fprintf('%-22s %10d %10d\n',    'N sections',       sum(slow),  sum(~slow));
fprintf('%-22s %10.1f %10.1f\n', 'median speed',     median(v(slow)),   median(v(~slow)));
fprintf('%-22s %10.0f %10.0f\n', 'median nFrames',   median(nFr(slow)), median(nFr(~slow)));
fprintf('%-22s %10.1f %10.1f\n', '% step reduced',   100*mean(red(slow)), 100*mean(red(~slow)));
fprintf('%-22s %10.1f %10.1f\n', '% on z>0 surface', 100*mean(sf(slow)==1), 100*mean(sf(~slow)==1));
fprintf('%-22s %10.3f %10.3f\n', 'median |k|',       median(k(slow)),   median(k(~slow)));

%% Refit on clean sections only
keep = nFr >= minFrames & ~red;
fprintf('\n--- Clean sections (nFrames >= %d, step not reduced): %d of %d ---\n', ...
    minFrames, sum(keep), numel(v));
kc = k(keep);  vc = v(keep);  tc = trk(keep);
pC = fitReport(vc, kc, tc, 'clean sections', nBoot);
for q = [1 5]
    m = vc >= prctile(vc, q);
    fitReport(vc(m), kc(m), tc(m), sprintf('clean, drop slowest %d%%', q), 0);
end

res = log(kc) - (pC(1) + pC(2)*log(vc));
fprintf('Scatter of log k around the robust fit: sd = %.2f (upper bound: includes estimator noise)\n', std(res));

%% ---------- Local functions ----------
function p = fitReport(v, k, trk, label, nBoot)
% Robust power-law fit (main) + OLS (comparison) + Spearman rho.
% nBoot > 0 adds a by-track bootstrap 95% CI for the robust exponent.
% Returns p = [intercept, slope] of the robust fit in log-log space.
lv = log(v);  lk = log(k);
p = robustfit(lv, lk);
pOLS = polyfit(lv, lk, 1);
rho = corr(v, k, 'Type', 'Spearman');
ciStr = '';
if nBoot > 0
    ci = bootRobustB(lv, lk, trk, nBoot);
    ciStr = sprintf('(95%% CI %+.2f to %+.2f)', ci);
end
fprintf('%-26s robust b = %+.2f %-24s | OLS b = %+.2f | rho = %+.3f | N = %d\n', ...
    label, p(2), ciStr, pOLS(1), rho, numel(v));
end

function ci = bootRobustB(lv, lk, trk, nBoot)
% Bootstrap over tracks (sections of one track are not independent).
[~, ~, g] = unique(trk);
members = accumarray(g, (1:numel(g))', [], @(x){x});
nT = numel(members);
bB = nan(nBoot,1);
for i = 1:nBoot
    j = vertcat(members{randi(nT, nT, 1)});
    pp = robustfit(lv(j), lk(j));
    bB(i) = pp(2);
end
ci = prctile(bB, [2.5 97.5]);
end
