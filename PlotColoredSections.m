clear; clc;

% --- Specify the results file ---
scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder), scriptFolder = pwd; end

[resultsFileName, resultsPath] = uigetfile(fullfile(scriptFolder, '*.mat'), ...
    'Select curvature results file');
if isequal(resultsFileName, 0)
    error('No results file selected. Aborting.');
end
resultsFile = fullfile(resultsPath, resultsFileName);

rawDataFile = 'Ref0722_unfiltered.mat';

smoothWindow = 1; % applied to k after loading
method = 'sgolay';

% --- Load results and raw trajectory data ---
resultsData = load(resultsFile);
allTrackResults = resultsData.allTrackResults;

rawData = load(rawDataFile);
S = rawData.AllBugs_unfilt;

% --- Specify which tracks to plot ---

trackIdsToPlot = [3461
1981
4057
4859
6803
4752
5846
1787
8472];
% rng(1)
% nTracksAvailable = numel(rawData.AllBugs_unfilt);
% trackIdsToPlot = randperm(nTracksAvailable, 10);

% --- Plot each requested track: smoothed curvature, with a trajectory inset ---
skippedTracks = [];

for i = 1:numel(trackIdsToPlot)
    trackIdx = trackIdsToPlot(i);

    try
        % --- Concatenate k, smoothed, with NaN gaps between sections ---
        sections = allTrackResults{trackIdx}.sections;
        if isempty(sections)
            error('No near-surface sections for this track.');
        end

        frameAll = [];
        kAll = [];
        for s = 1:numel(sections)
            kSmoothed = smoothCurvature(sections(s).k, smoothWindow, method);
            frameAll = [frameAll; sections(s).frame(:); NaN]; %#ok<AGROW>
            kAll = [kAll; kSmoothed(:); NaN]; %#ok<AGROW>
        end
        if ~isempty(frameAll)
            frameAll(end) = []; kAll(end) = [];
        end

         % --- Pull the raw trajectory for the inset ---
        M = rawData.AllBugs_unfilt{trackIdx};
        validRows = ~any(isnan(M(:,2:4)), 2);
        M = M(validRows, :);
        xTraj = M(:,2);
        yTraj = M(:,3);
        frameTraj = M(:,1);

        % --- Main curvature plot ---
        figure;
        mainAx = axes;
        cols = lines(numel(sections));   % one color per section, reused in the inset
        hold(mainAx, 'on');
        nRed = 0;
        for s = 1:numel(sections)
            kS = smoothCurvature(sections(s).k, smoothWindow, method);
            if isfield(sections, 'stepReduced') && sections(s).stepReduced
                lineStyle = '--';  nRed = nRed + 1;
            else
                lineStyle = '-';
            end
            plot(mainAx, sections(s).frame, kS, lineStyle, 'Color', cols(s,:), 'LineWidth', 1.5);
        end
        hold(mainAx, 'off');
        if nRed > 0
            text(mainAx, 0.02, 0.95, sprintf('step reduced to %d in %d/%d sections (dashed)', ...
                sections(find([sections.stepReduced],1)).stepUsed, nRed, numel(sections)), ...
                'Units', 'normalized', 'Color', 'k', 'VerticalAlignment', 'top');
        end
       
        maxAbsK = max(abs(kAll), [], 'omitnan');
        yBound = max(0.3, maxAbsK);
        ylim(mainAx, [-yBound, yBound]);

        xlabel(mainAx, 'frame');
        ylabel(mainAx, 'curvature');
        title(mainAx, sprintf('Track %d: smoothed curvature (window=%d)', trackIdx, smoothWindow));
        grid(mainAx, 'on');

        % --- Trajectory inset, showing ONLY the near-surface sections
        % that curvature was actually computed on, colored by frame
        % (time) so you can see the direction/progression of motion.
        % Color scale is fixed to the FULL track's frame range, so
        % absolute time position is comparable across sections, not
        % just progression within each section individually. ---
        insetAx = axes('Position', [0.68 0.68 0.2 0.2]);
        hold(insetAx, 'on');
        for s = 1:numel(sections)
            idxSec = sections(s).startRow : sections(s).endRow;
            plot(insetAx, xTraj(idxSec), yTraj(idxSec), '.-', 'Color', cols(s,:), ...
                'MarkerSize', 8, 'LineWidth', 0.5);
        end
        hold(insetAx, 'off');
        axis(insetAx, 'equal');
        set(insetAx, 'XTick', [], 'YTick', [], 'Box', 'on');
        

    catch ME
        skippedTracks(end+1,1) = trackIdx; %#ok<AGROW>
        fprintf('Skipped track %d: %s\n', trackIdx, ME.message);
    end
end