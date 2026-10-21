function [multBest, multGrid, rmse, bias] = sweep_tolMult(trajNoisy, kTrue, methodFcn, tFrames, multGrid, useSignedComparison)
% methodFcn(traj, t, tolMult) -> [kEst, idx]
% trajNoisy : cell array, each cell is N x 2 [x y] noisy track
% kTrue     : cell array, true curvature per track (N x 1)
% tFrames   : cell array, noise-free parameter per track (frame numbers)
% multGrid  : candidate tolerance multipliers (optional)
% useSignedComparison : if true, compares signed kEst vs signed kTrue
%                        (only meaningful for shapes with signed kTrue,
%                        e.g. SpatiallyVarying). Default false.

if nargin < 4, tFrames = cell(size(trajNoisy)); end
if nargin < 5 || isempty(multGrid), multGrid = logspace(-1, 1, 20); end
if nargin < 6, useSignedComparison = false; end

rmse = nan(size(multGrid));
bias = nan(size(multGrid));

for j = 1:numel(multGrid)
    err = [];
    for k = 1:numel(trajNoisy)
        traj = trajNoisy{k};
        if size(traj,1) < 5, continue; end

        t = [];
        if k <= numel(tFrames) && ~isempty(tFrames{k})
            t = tFrames{k};
        end

        [kEst, idx] = methodFcn(traj, t, multGrid(j));
        if isempty(kEst), continue; end

        kt = kTrue{k}(idx(:));
        if useSignedComparison
            err = [err; kEst(:) - kt]; %#ok<AGROW>
        else
            err = [err; abs(kEst(:)) - abs(kt)]; %#ok<AGROW>
        end
    end
    rmse(j) = sqrt(mean(err.^2, 'omitnan'));
    bias(j) = mean(err, 'omitnan');
end

[~, iBest] = min(rmse);
multBest = multGrid(iBest);

figure;
subplot(2,1,1); semilogx(multGrid, rmse, '-o'); hold on;
semilogx(multBest, rmse(iBest), 'r*', 'MarkerSize', 12);
xlabel('tolerance multiplier'); ylabel('RMSE'); grid on;
subplot(2,1,2); semilogx(multGrid, bias, '-o'); yline(0, 'k--');
xlabel('tolerance multiplier'); ylabel('Bias (est - true)'); grid on;

fprintf('Best multiplier = %.4g (RMSE = %.4g)\n', multBest, rmse(iBest));
end