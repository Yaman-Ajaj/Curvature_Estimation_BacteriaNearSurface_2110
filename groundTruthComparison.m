clear; 
clc;

scriptFolder = fileparts(mfilename('fullpath'));
if isempty(scriptFolder), scriptFolder = pwd; end

[fileNames, filePath] = uigetfile(fullfile(scriptFolder, '*.mat'), ...
    'Select result files (Ctrl/Shift-click for multiple)', ...
    'MultiSelect', 'on');

if isequal(fileNames, 0)
    error('No files selected. Aborting.');
end

% With MultiSelect, a single selection returns a char instead of a cell
if ischar(fileNames)
    fileNames = {fileNames};
end

resultsFiles = fullfile(filePath, fileNames(:)); % column cell array of full paths


syntheticFile = 'Fake_synthetic_trajectories_Clean.mat';
syntheticData = load(syntheticFile);
trajectoryType = syntheticData.trajectoryType;
trajectoryRadii = syntheticData.trajectoryRadii;
trajectorySemiAxes = syntheticData.trajectorySemiAxes;
trajectoryTheta = syntheticData.trajectoryTheta; % requires regenerated synthetic file

allAbsError = [];
allRelError = [];
allFileLabels = [];
allShapeLabels = {};

for f = 1:numel(resultsFiles)
    thisFile = resultsFiles{f};
    if ~isfile(thisFile)
        fprintf('File not found, skipping: %s\n', thisFile);
        continue;
    end

    resultsData = load(thisFile);
    if ~isfield(resultsData, 'allTrackResults')
        fprintf('File missing allTrackResults, skipping: %s\n', thisFile);
        continue;
    end
    allTrackResults = resultsData.allTrackResults;

    absErrThisFile = [];
    relErrThisFile = [];
    shapeThisFile = {};

    for t = 1:numel(allTrackResults)
        if isempty(allTrackResults{t})
            continue;
        end

        thisType = trajectoryType{t};
        sections = allTrackResults{t}.sections;

        for s = 1:numel(sections)
            kVals = sections(s).k;
            idxVals = sections(s).idx; % row indices into the ORIGINAL trajectory
            validMask = isfinite(kVals);
            kVals = kVals(validMask);
            idxVals = idxVals(validMask);

            if isempty(kVals)
                continue;
            end

            % --- Trim a percentage of points at each end of the RUN ---
            runLen = sections(s).endRow - sections(s).startRow + 1;
            pct = edgeTrimPercent(runLen);
            nTrim = round(pct/100 * runLen);

            if nTrim > 0
                keepMask = (idxVals > sections(s).startRow + nTrim - 1) & ...
                   (idxVals < sections(s).endRow   - nTrim + 1);
                kVals   = kVals(keepMask);
                idxVals = idxVals(keepMask);
            end

            if isempty(kVals)
                continue;
            end

            % --- Compute the TRUE curvature per point, depending on shape ---
            switch thisType
                case {'Circular', 'GentleArc'}
                    % Constant true curvature everywhere: 1/radius
                    trueK = (1 / trajectoryRadii(t)) * ones(size(kVals));

                case 'Elliptical'
                    % Position-dependent true curvature via the exact
                    % ellipse curvature formula, using each point's
                    % actual theta (saved by the generator, not
                    % re-derived) — kappa(theta) = a*b / (a^2*sin^2(theta) + b^2*cos^2(theta))^(3/2)
                    a = trajectorySemiAxes(t,1);
                    b = trajectorySemiAxes(t,2);
                    thetaAll = trajectoryTheta{t};
                    thetaVals = thetaAll(idxVals);
                    trueK = (a*b) ./ (a^2*sin(thetaVals).^2 + b^2*cos(thetaVals).^2).^1.5;
                    trueK = trueK(:);

                otherwise
                    continue; % unknown type, skip rather than guess
            end

            absErr = abs(kVals(:) - trueK(:));
            relErr = 100 * absErr ./ trueK(:);

            absErrThisFile = [absErrThisFile; absErr]; %#ok<AGROW>
            relErrThisFile = [relErrThisFile; relErr]; %#ok<AGROW>
            shapeThisFile = [shapeThisFile; repmat({thisType}, numel(absErr), 1)]; %#ok<AGROW>
        end
    end

    if isempty(absErrThisFile)
        fprintf('File has no usable sections, skipping: %s\n', thisFile);
        continue;
    end

    allAbsError = [allAbsError; absErrThisFile]; %#ok<AGROW>
    allRelError = [allRelError; relErrThisFile]; %#ok<AGROW>
    allFileLabels = [allFileLabels; repmat(f, numel(absErrThisFile), 1)]; %#ok<AGROW>
    allShapeLabels = [allShapeLabels; shapeThisFile]; %#ok<AGROW>
end

% --- Print summary stats per file, split by shape type ---
uniqueFiles = unique(allFileLabels);
uniqueShapes = unique(allShapeLabels);

fprintf('--- Relative error (%%) by file and shape type ---\n');
fprintf('%-45s %-12s %8s %10s %10s %10s\n', 'resultsFile', 'shape', 'N', 'Median', 'Mean', 'p90');
for fi = 1:numel(uniqueFiles)
    for si = 1:numel(uniqueShapes)
        mask = (allFileLabels == uniqueFiles(fi)) & strcmp(allShapeLabels, uniqueShapes{si});
        vals = allRelError(mask);
        if isempty(vals)
            continue;
        end
        fprintf('%-45s %-12s %8d %10.2f %10.2f %10.2f\n', ...
            resultsFiles{uniqueFiles(fi)}, uniqueShapes{si}, numel(vals), median(vals), mean(vals), prctile(vals,90));
    end
end