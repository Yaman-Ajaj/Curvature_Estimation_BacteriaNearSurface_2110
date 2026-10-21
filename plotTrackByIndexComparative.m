clear; clc;
% --- Specify the two saved results files to compare ---
groundTruthFile = 'groundTruth.mat';
% --- Pick the result file via a popup ---
scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder)
    scriptFolder = pwd; % fallback when running as a section/selection rather than a saved script
end
[fileName, filePath] = uigetfile(fullfile(scriptFolder, '*.mat'), ...
    'Select the curvature results file to compare');
if isequal(fileName, 0)
    error('No file selected. Aborting.');
end
resultFile = fullfile(filePath, fileName);
fprintf('Comparing against: %s\n', fileName);

% --- Specify which tracks to compare ---
%trackIdsToPlot = [8, 27, 9, 24, 17, 18, 4, 5, 28, 30];
rng(3)
trackIdsToPlot = randperm(30, 5);
smoothWindow = 1;        % applied to both methods' k after loading, for a fair comparison
minYLimits = [-0.4 0.4]; % default range; expands only if a track exceeds it
useAbsForTangent = false; % set true to compare magnitudes (then negative half of axis will be empty for that curve)

% --- Load both results files ---
groundTruthData = load(groundTruthFile);
allTrackResults = groundTruthData.allTrackResults;
resultsData = load(resultFile);
allTrackResultsAnalyzed = resultsData.allTrackResults;

% --- Plot each requested track: circle-fit vs tangent-angle, both smoothed ---
skippedTracks = [];
for i = 1:numel(trackIdsToPlot)
    trackIdx = trackIdsToPlot(i);
    try
        % --- Concatenate circle-fit k, smoothed, with NaN gaps between sections ---
        circleSections = allTrackResults{trackIdx}.sections;
        if isempty(circleSections)
            error('No near-surface sections for this track (circle-fit results).');
        end
        frameCircleAll = [];
        kCircleAll = [];
        for s = 1:numel(circleSections)
            kSmoothed = smoothCurvature(circleSections(s).k, smoothWindow);
            frameCircleAll = [frameCircleAll; circleSections(s).frame(:); NaN]; %#ok<AGROW>
            kCircleAll = [kCircleAll; kSmoothed(:); NaN]; %#ok<AGROW>
        end
        if ~isempty(frameCircleAll)
            frameCircleAll(end) = []; kCircleAll(end) = [];
        end

        % --- Concatenate tangent-angle k, smoothed, with NaN gaps between sections ---
        tangentSections = allTrackResultsAnalyzed{trackIdx}.sections;
        if isempty(tangentSections)
            error('No near-surface sections for this track (tangent-angle results).');
        end
        frameTangentAll = [];
        kTangentAll = [];
        for s = 1:numel(tangentSections)
            kRaw = tangentSections(s).k;
            if useAbsForTangent
                kRaw = abs(kRaw);
            end
            kSmoothed = smoothCurvature(kRaw, smoothWindow);
            frameTangentAll = [frameTangentAll; tangentSections(s).frame(:); NaN]; %#ok<AGROW>
            kTangentAll = [kTangentAll; kSmoothed(:); NaN]; %#ok<AGROW>
        end
        if ~isempty(frameTangentAll)
            frameTangentAll(end) = []; kTangentAll(end) = [];
        end

        % --- Plot both on the same axes ---
        figure;
        plot(frameCircleAll, kCircleAll, 'b-', 'LineWidth', 1.2, 'DisplayName', 'Ground Truth');
        hold on;
        plot(frameTangentAll, kTangentAll, 'r-', 'LineWidth', 1.2, 'DisplayName', 'Analyzed Curvature');
        yline(0, 'k:', 'HandleVisibility', 'off'); % zero reference line
        hold off;

        % --- Y-axis: default +/-0.4, expanded only if this track exceeds it ---
        allK = [kCircleAll; kTangentAll];
        yMin = min(minYLimits(1), min(allK, [], 'omitnan'));
        yMax = max(minYLimits(2), max(allK, [], 'omitnan'));
        ylim([yMin yMax]);

        xlabel('frame');
        ylabel('curvature');
        title(sprintf('Track %d: curvature method comparison', trackIdx));
        legend('show');
        grid on;

    catch ME
        skippedTracks(end+1,1) = trackIdx; %#ok<AGROW>
        fprintf('Skipped track %d: %s\n', trackIdx, ME.message);
    end
end