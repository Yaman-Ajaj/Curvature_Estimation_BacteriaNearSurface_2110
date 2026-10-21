clear; 
%clc; 
%close all

% --- Load data ---
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

% --- Choose method ---
methodOptions = {
    'Spline-smoothed + circle-fit'
    'Spline-smoothed + least-squares circle'
    'Spline-smoothed + tangent-angle (signed)'
    'Spline-smoothed derivative (signed)'
};
[mSel, ok] = listdlg('PromptString','Select method:', ...
    'SelectionMode','single','ListString',methodOptions,'ListSize',[300,120]);
if ~ok, error('No method selected.'); end

% --- Ask for parameters ---
usesStep = (mSel ~= 4);   % spline-derivative has no step
if usesStep
    prompt = {'sigma (noise std):', 'tolMult:', 'step / windowSize:'};
    defaults = {'0.2', '1.29', '1'};
else
    prompt = {'sigma (noise std):', 'tolMult:'};
    defaults = {'0.2', '1.29'};
end
answer = inputdlg(prompt, 'Method parameters', [1 40], defaults);
if isempty(answer), error('Cancelled.'); end

sigma   = str2double(answer{1});
tolMult = str2double(answer{2});
if usesStep, step = str2double(answer{3}); end

% --- Build the method handle ---
switch mSel
    case 1
        methodFcn = @(traj,t) curvature_estimation_splineSmoothed_circlefit(traj, sigma, step, t, tolMult);
    case 2
        methodFcn = @(traj,t) curvature_estimation_splineSmoothed_leastSquaresCircle(traj, sigma, step, t, tolMult);
    case 3
        methodFcn = @(traj,t) curvature_estimation_splineSmoothed(traj, sigma, step, t, tolMult);
    case 4
        methodFcn = @(traj,t) curvature_estimation_splineSmoothed_derivative(traj, sigma, t, tolMult);
end

% --- Signed comparison only makes sense for the signed method on signed shapes ---
signedShapes = {'SpatiallyVarying','FastVarying'};
useSigned = (mSel == 4) && ismember(restrictedType, signedShapes);

% --- Evaluate ---
err = [];
for k = 1:numel(trajNoisy)
    traj = trajNoisy{k};
    if size(traj,1) < 5, continue; end

    t = [];
    if k <= numel(tFrames) && ~isempty(tFrames{k}), t = tFrames{k}; end

    [kEst, idx] = methodFcn(traj, t);
    if isempty(kEst), continue; end

    kt = kTrue{k}(idx(:));
    if useSigned
        err = [err; kEst(:) - kt]; %#ok<AGROW>
    else
        err = [err; abs(kEst(:)) - abs(kt)]; %#ok<AGROW>
    end
end

rmse    = sqrt(mean(err.^2, 'omitnan'));
bias    = mean(err, 'omitnan');
nPoints = sum(~isnan(err));

% --- Report ---
fprintf('\n--- Result ---\n');
fprintf('Method     : %s\n', methodOptions{mSel});
fprintf('Shape type : %s\n', restrictedType);
if usesStep
    fprintf('Parameters : sigma = %g, tolMult = %g, step = %g\n', sigma, tolMult, step);
else
    fprintf('Parameters : sigma = %g, tolMult = %g\n', sigma, tolMult);
end
fprintf('Comparison : %s\n', ternary(useSigned, 'signed', 'unsigned (magnitude)'));
fprintf('RMSE       : %.5f\n', rmse);
fprintf('Bias       : %+.5f\n', bias);
fprintf('N points   : %d\n', nPoints);

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end