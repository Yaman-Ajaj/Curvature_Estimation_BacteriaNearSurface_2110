clear; clc;

%% Files
filteredFile   = 'Ref0722_filtered.mat';     % used ONLY to detect near-surface runs
unfilteredFile = 'Ref0722_unfiltered.mat';   % <-- put your unfiltered file name here

fData = load(filteredFile);
uData = load(unfilteredFile);
SF = fData.AllBugs_filt;      % adjust if the variable layout differs
SU = uData.AllBugs_unfilt;

nTracks = min(numel(SF), numel(SU));
if numel(SF) ~= numel(SU)
    warning('Track counts differ (filtered %d, unfiltered %d). Using the first %d.', ...
        numel(SF), numel(SU), nTracks);
end

%% Z-boundary parameters (applied to the FILTERED z)
surfaceZmin = 35.9535;
surfaceZmax = 35.8652;
tolerance   = 3.5;

%% Method selection
methodOptions = {
    'Spline-smoothed + circle-fit'
    'Spline-smoothed + least-squares circle'
    'Spline-smoothed + tangent-angle (signed)'
    'Spline-smoothed derivative (signed)'
    'Circle-fit (osculating circle)'
};
[selection, ok] = listdlg('PromptString', 'Select curvature method:', ...
    'SelectionMode', 'single', 'ListString', methodOptions, 'ListSize', [300, 120]);
if ~ok, error('No method selected. Aborting.'); end

%% Parameters (update from the new sweeps)
sigma         = 0.25;   % effective noise of the unfiltered data
tolMult       = 1.29;   % <-- from new sweep
curvatureStep = 7;      % <-- from new sweep (circle / tangent-angle)
lsWindowSize  = 9;      % <-- from new sweep (least-squares circle)
%minRunLength  = 2*curvatureStep + 3;   % shorter runs cannot give curvature

sg = strrep(num2str(sigma), '.', '');
tm = strrep(num2str(tolMult), '.', '');

switch selection
    case 1
        curvatureMethod = sprintf('Real_splinesmoothed_circlefit_sig%s_m%s', sg, tm);
        curvatureFcn = @(x,y,z,t,p) curvature_estimation_splineSmoothed_circlefit([x, y], sigma, p, t, tolMult);
        paramFull = curvatureStep;  paramFallback = 1;  minLenFcn = @(p) 2*p + 3;
    case 2
        curvatureMethod = sprintf('Real_splinesmoothed_leastSquaresCircle_sig%s_m%s_w%d', sg, tm, lsWindowSize);
        curvatureFcn = @(x,y,z,t,p) curvature_estimation_splineSmoothed_leastSquaresCircle([x, y], sigma, p, t, tolMult);
        paramFull = lsWindowSize;  paramFallback = 5;  minLenFcn = @(p) p + 2;
    case 3
        curvatureMethod = sprintf('Real_splinesmoothed_tangentAngle_sig%s_m%s', sg, tm);
        curvatureFcn = @(x,y,z,t,p) curvature_estimation_splineSmoothed([x, y], sigma, p, t, tolMult);
        paramFull = curvatureStep;  paramFallback = 1;  minLenFcn = @(p) 2*p + 3;
    case 4
        curvatureMethod = sprintf('Real_splineDerivative_sig%s_m%s', sg, tm);
        curvatureFcn = @(x,y,z,t,p) curvature_estimation_splineSmoothed_derivative([x, y], sigma, t, tolMult);
        paramFull = NaN;  paramFallback = NaN;  minLenFcn = @(p) 5;   % no step
    case 5
        curvatureMethod = 'Real_circlefit';
        curvatureFcn = @(x,y,z,t,p) deal_kidx([x,y], p);
        paramFull = curvatureStep;  paramFallback = 1;  minLenFcn = @(p) 2*p + 3;
end
fprintf('Using method: %s\n', curvatureMethod);

%% Main loop
allTrackResults = cell(nTracks, 1);
minOut = 5;   % define once, above the track loop

for k = 1:nTracks

    MF = SF{k};   MU = SU{k};
    if isempty(MF) || isempty(MU) || size(MF,2) < 4 || size(MU,2) < 4
        continue;
    end

    % ---- Filtered data: detect near-surface runs ----
    MF = MF(~any(isnan(MF(:,2:4)), 2), :);
    if isempty(MF), continue; end

    frameF = MF(:,1);
    zF     = MF(:,4);
    nearSurface = (abs(zF - surfaceZmax) <= tolerance) | ...
                  (abs(zF + surfaceZmin) <= tolerance);

    % split on surface departures AND on frame gaps
    segId  = cumsum([1; diff(frameF) > 1]);
    runKey = double(nearSurface) .* segId;
    runStarts = find([true; diff(runKey) ~= 0]);
    runEnds   = [runStarts(2:end) - 1; numel(runKey)];
    isRunNearSurface = runKey(runStarts) > 0;

    % ---- Unfiltered data: valid rows only; sections index into THIS matrix ----
    MU = MU(~any(isnan(MU(:,2:4)), 2), :);
    if isempty(MU), continue; end
    frameU = MU(:,1);

    trackSections = struct([]);

    for r = 1:numel(runStarts)
        if ~isRunNearSurface(r), continue; end

        % frames of this run according to the filtered data
        runFrames = frameF(runStarts(r):runEnds(r));

        % same frames in the unfiltered data
        rowsU = find(ismember(frameU, runFrames));
        if numel(rowsU) < minLenFcn(paramFallback), continue; end

        % unfiltered rows missing a frame -> split so no calculation bridges a gap
        subSeg = cumsum([1; diff(frameU(rowsU)) > 1]);

        for sg_i = unique(subSeg)'
            idxRun = rowsU(subSeg == sg_i);

            paramUsed = paramFull;
            stepReduced = false;
            if numel(idxRun) < minLenFcn(paramFull)
                
                paramUsed = max(1, min(paramFull, floor((numel(idxRun) - minOut)/2)));
                stepReduced = paramUsed < paramFull;
                if numel(idxRun) < minLenFcn(paramFallback), continue; end   % too short even for the fallback
            end

            xRun = MU(idxRun,2);
            yRun = MU(idxRun,3);
            zRun = MU(idxRun,4);
            frameRun = frameU(idxRun);

            [kRun, idxLocal] = curvatureFcn(xRun, yRun, zRun, frameRun, paramUsed);
            if isempty(kRun), continue; end

            kFin = kRun(isfinite(kRun));

            secResult = struct();
            secResult.sectionIndex = numel(trackSections) + 1;
            secResult.startRow   = idxRun(1);     % rows of the NaN-filtered UNFILTERED matrix
            secResult.endRow     = idxRun(end);
            secResult.startFrame = frameRun(1);
            secResult.endFrame   = frameRun(end);
            secResult.frame      = frameRun(idxLocal);
            secResult.k          = kRun;
            secResult.idx        = idxLocal;      % local to this run
            secResult.meanCurvature = mean(kFin);
            secResult.stdCurvature  = std(kFin);
            secResult.stepUsed    = paramUsed;
            secResult.stepReduced = stepReduced;

            trackSections = [trackSections; secResult]; %#ok<AGROW>
        end
    end

    allTrackResults{k} = struct('trackIndex', k, 'sections', trackSections);
end

%% Save
resultsFile = sprintf('curvatureResults_%s_step%d_%s.mat', curvatureMethod, ...
    curvatureStep, char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
save(resultsFile, 'allTrackResults', 'surfaceZmin', 'surfaceZmax', 'tolerance', ...
    'curvatureStep', 'curvatureMethod', 'sigma', 'tolMult', 'filteredFile', 'unfilteredFile');
fprintf('Results saved to %s\n', resultsFile);