function [pBest, pGrid, rmse, bias] = sweep_splineP(trajNoisy, kTrue, step, pGrid, tFrames)
% trajNoisy : cell array, each cell is N x 2 [x y] noisy track
% kTrue     : cell array, true curvature at each point of each track (N x 1)
% step      : curvature step (e.g. 1)
% pGrid     : candidate p values (optional)
% tFrames   : cell array, noise-free parameter per track (e.g. frame
%             numbers), same length as trajNoisy{k}. Optional; if omitted
%             or a cell is empty, falls back to cumulative chord length
%             for that track.

if nargin < 4 || isempty(pGrid)
    pGrid = logspace(-4, 0, 30);
end
if nargin < 5
    tFrames = cell(size(trajNoisy));
end

rmse = nan(size(pGrid));
bias = nan(size(pGrid));

for j = 1:numel(pGrid)
    err = [];
    for k = 1:numel(trajNoisy)
        traj = trajNoisy{k};
        if size(traj,1) < 5, continue; end

        t = [];
        if k <= numel(tFrames) && ~isempty(tFrames{k})
            t = tFrames{k};
        end

        [kEst, idx] = curvature_estimation_splineSmoothed_circlefit(traj, pGrid(j), step, t);
        err = [err; kEst(:) - abs(kTrue{k}(idx(:)))]; %#ok<AGROW>
    end
    rmse(j) = sqrt(mean(err.^2, 'omitnan'));
    bias(j) = mean(err, 'omitnan');
end

[~, iBest] = min(rmse);
pBest = pGrid(iBest);

figure;
subplot(2,1,1); semilogx(pGrid, rmse, '-o'); hold on;
semilogx(pBest, rmse(iBest), 'r*', 'MarkerSize', 12);
xlabel('p'); ylabel('RMSE'); grid on;
subplot(2,1,2); semilogx(pGrid, bias, '-o'); yline(0, 'k--');
xlabel('p'); ylabel('Bias (est - true)'); grid on;

fprintf('Best p = %.4g (RMSE = %.4g)\n', pBest, rmse(iBest));
end