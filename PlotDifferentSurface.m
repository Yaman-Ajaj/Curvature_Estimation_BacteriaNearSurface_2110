clear; clc;

%% Settings
scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder), scriptFolder = pwd; end
[fn, fp] = uigetfile(fullfile(scriptFolder,'*.mat'), 'Select curvature results file');
if isequal(fn,0), error('No file selected.'); end
R = load(fullfile(fp,fn));
allTrackResults = R.allTrackResults;
SF = load(R.filteredFile).AllBugs_filt;   % same filtered file the analysis used
SU = load(R.unfilteredFile).AllBugs_unfilt;

smoothWindow = 1;  smoothMethod = 'sgolay';
colPos = [0.85 0.33 0.10];    % z > 0 surface
colNeg = [0.00 0.45 0.74];    % z < 0 surface
labels = {'z > 0 surface', 'z < 0 surface'};

% rng(2);
% trackIds = randperm(numel(allTrackResults), 10);
trackIds = [136,146,149,168,175,189,201,218,229];

for t = trackIds
    if isempty(allTrackResults{t}) || isempty(allTrackResults{t}.sections), continue; end
    MF = SF{t};  MF = MF(~any(isnan(MF(:,2:4)),2), :);
    MU = SU{t};  MU = MU(~any(isnan(MU(:,2:4)),2), :);
    secs = allTrackResults{t}.sections;

    figure; hold on;
    mainAx = gca;
    maxAbsK = 0;
    hPos = plot(NaN,NaN,'-','Color',colPos,'LineWidth',1.5);   % legend handles
    hNeg = plot(NaN,NaN,'-','Color',colNeg,'LineWidth',1.5);
    secCol = zeros(numel(secs),3);
    for s = 1:numel(secs)
        inSec = MF(:,1) >= secs(s).startFrame & MF(:,1) <= secs(s).endFrame;
        c = colPos;
        if median(MF(inSec,4)) < 0, c = colNeg; end
        secCol(s,:) = c;
        kS = smoothCurvature(secs(s).k, smoothWindow, smoothMethod);
        plot(mainAx, secs(s).frame, kS, '-', 'Color', c, 'LineWidth', 1.5);
        maxAbsK = max(maxAbsK, max(abs(kS), [], 'omitnan'));
    end

    yBound = max(0.3, maxAbsK);
    ylim(mainAx, [-yBound, yBound]);

    yline(0,'k:','HandleVisibility','off');
    xlabel('frame'); ylabel('curvature'); grid on;
    title(sprintf('Track %d: curvature by surface', t));
    legend([hPos hNeg], labels, 'Location', 'northwest');

    insetAx = axes('Position', [0.68 0.68 0.2 0.2]);
    hold(insetAx, 'on');
    for s = 1:numel(secs)
        idxSec = secs(s).startRow : secs(s).endRow;
        plot(insetAx, MU(idxSec,2), MU(idxSec,3), '.-', 'Color', secCol(s,:), ...
            'MarkerSize', 8, 'LineWidth', 0.5);
    end
    hold(insetAx, 'off');
    axis(insetAx, 'equal');
    set(insetAx, 'XTick', [], 'YTick', [], 'Box', 'on');
end