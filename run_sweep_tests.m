clear; 
%clc; 
%close all
% --- Load data (always needed) ---
D = load('Fake_synthetic_trajectories_Noise.mat');
kTrue = getTrueCurvature('Fake_synthetic_trajectories_Clean.mat');
trajNoisy = cellfun(@(M) M(:,2:3), D.trajectories, 'UniformOutput', false);
tFrames   = cellfun(@(M) M(:,1),   D.trajectories, 'UniformOutput', false);

isClosed = ismember(D.trajectoryType, {'Circular','Elliptical'});
for k = find(isClosed)'
    trajNoisy{k}(end,:) = [];
    kTrue{k}(end)       = [];
    tFrames{k}(end)     = [];
end

sigma = 0.24;

% --- Restrict to a shape type (optional) ---
excludedTypes = {'SpatiallyVarying','FastVarying'};
typeOptions = [{'All'; 'All (excluding SpatiallyVarying & FastVarying)'}; unique(D.trajectoryType)];
[sel, ok] = listdlg('PromptString','Restrict to shape type:', ...
    'SelectionMode','single','ListString',typeOptions,'ListSize',[350,160]);
if ~ok, error('No selection made.'); end
restrictedType = typeOptions{sel};

switch sel
    case 1   % All
        mask = true(size(D.trajectoryType));
    case 2   % All except the two varying-curvature types
        mask = ~ismember(D.trajectoryType, excludedTypes);
    otherwise
        mask = strcmp(D.trajectoryType, restrictedType);
end

trajNoisy = trajNoisy(mask);
kTrue     = kTrue(mask);
tFrames   = tFrames(mask);


% Use signed comparison only when the restricted type has signed kTrue
% (currently just SpatiallyVarying). Circles/ellipses/arcs are unsigned
% by construction, so signed comparison there would be meaningless.
useSigned = strcmp(restrictedType, 'SpatiallyVarying');

% --- Choose which sweep test to run ---
testOptions = {
    'Spline tolerance (spline-derivative method)'
    'Spline tolerance x step (circle-fit method)'
    'Spline tolerance x step (least-squares circle method)'
    'Spline tolerance x step (tangent-angle method)'
};
[testSel, ok] = listdlg('PromptString','Select sweep test:', ...
    'SelectionMode','single','ListString',testOptions,'ListSize',[350,120]);
if ~ok, error('No test selected.'); end

switch testSel
    case 1
        methodFcn = @(traj,t,tolMult) curvature_estimation_splineSmoothed_derivative(traj, sigma, t, tolMult);
        sweep_tolMult(trajNoisy, kTrue, methodFcn, tFrames, [], useSigned);

    case 2
        stepGrid = [1 2 3 5 7 10 15 20 25 30 35];
        methodFcn = @(traj,t,tolMult,step) curvature_estimation_splineSmoothed_circlefit(traj, sigma, step, t, tolMult);
        sweep_splineCirclefit_2D(trajNoisy, kTrue, tFrames, methodFcn, [], stepGrid, false);
        

    case 3
        stepGridOdd = [5 7 9 11 15 21 31];  % all odd, matches windowSize requirement directly
        methodFcn = @(traj,t,tolMult,step) curvature_estimation_splineSmoothed_leastSquaresCircle(traj, sigma, step, t, tolMult);
        sweep_splineCirclefit_2D(trajNoisy, kTrue, tFrames, methodFcn, [], stepGridOdd, false);
    
    case 4
        stepGrid = [1 2 3 5 7 10 15 20 25 30 35];
        methodFcn = @(traj,t,tolMult,step) curvature_estimation_splineSmoothed(traj, sigma, step, t, tolMult);
        sweep_splineCirclefit_2D(trajNoisy, kTrue, tFrames, methodFcn, [], stepGrid, false);
end
