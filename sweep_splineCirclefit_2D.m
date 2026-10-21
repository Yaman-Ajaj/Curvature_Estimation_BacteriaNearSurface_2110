function [tolBest, stepBest, RMSE, tolGrid, stepGrid] = sweep_splineCirclefit_2D(trajNoisy, kTrue, tFrames, methodFcn, tolGrid, stepGrid, useSignedComparison)
% methodFcn(traj, t, tolMult, step) -> [kEst, idx]

if nargin < 5 || isempty(tolGrid),  tolGrid  = logspace(-1, 1, 10); end
if nargin < 6 || isempty(stepGrid), stepGrid = [1 2 3 5 7 10];      end
if nargin < 7, useSignedComparison = false; end

nT = numel(tolGrid); nS = numel(stepGrid);
RMSE = nan(nT, nS);

for it = 1:nT
    for is = 1:nS
        err = [];
        for k = 1:numel(trajNoisy)
            traj = trajNoisy{k};
            if size(traj,1) < 2*stepGrid(is)+3, continue; end
            tk = [];
            if k <= numel(tFrames) && ~isempty(tFrames{k}), tk = tFrames{k}; end

            [kEst, idx] = methodFcn(traj, tk, tolGrid(it), stepGrid(is));
            if isempty(kEst), continue; end

            kt = kTrue{k}(idx(:));
            if useSignedComparison
                err = [err; kEst(:) - kt]; %#ok<AGROW>
            else
                err = [err; abs(kEst(:)) - abs(kt)]; %#ok<AGROW>
            end
        end
        RMSE(it,is) = sqrt(mean(err.^2, 'omitnan'));
    end
end

[minPerCol, rowIdx] = min(RMSE, [], 1, 'omitnan');
[~, colIdx] = min(minPerCol);
tolBest = tolGrid(rowIdx(colIdx));
stepBest = stepGrid(colIdx);

figure('Color','w');
imagesc(stepGrid, log10(tolGrid), RMSE);
set(gca, 'YDir', 'normal');
xlabel('step'); ylabel('log_{10}(tolMult)');
title('Curvature RMSE');
colorbar;
xticks(stepGrid);
hold on;
plot(stepBest, log10(tolBest), 'r*', 'MarkerSize', 14, 'LineWidth', 2);

fprintf('Best: tolMult = %.3g, step = %d (RMSE = %.4g)\n', tolBest, stepBest, RMSE(rowIdx(colIdx), colIdx));
end