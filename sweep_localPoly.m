function [wBest, oBest, RMSE, windowSizes, polyOrders] = sweep_localPoly(trajNoisy, kTrue, tFrames, methodFcn, windowSizes, polyOrders, useSignedComparison)
% methodFcn(traj, t, windowSize, polyOrder) -> [kEst, idx]

if nargin < 4 || isempty(methodFcn)
    methodFcn = @(traj,t,w,o) curvature_estimation_localPoly(traj, w, o, t);
end
if nargin < 5 || isempty(windowSizes), windowSizes = [5 7 9 11 15 21 31]; end
if nargin < 6 || isempty(polyOrders),  polyOrders  = 2:6; end
if nargin < 7, useSignedComparison = false; end

nW = numel(windowSizes);
nP = numel(polyOrders);
RMSE = nan(nW, nP);

for iw = 1:nW
    w = windowSizes(iw);
    for ip = 1:nP
        o = polyOrders(ip);
        if o >= w
            continue;
        end
        err = [];
        for k = 1:numel(trajNoisy)
            traj = trajNoisy{k};
            if size(traj,1) < w, continue; end

            tk = [];
            if k <= numel(tFrames) && ~isempty(tFrames{k})
                tk = tFrames{k};
            end

            [kEst, idx] = methodFcn(traj, tk, w, o);
            if isempty(kEst), continue; end

            kt = kTrue{k}(idx(:));
            if useSignedComparison
                err = [err; kEst(:) - kt]; %#ok<AGROW>
            else
                err = [err; abs(kEst(:)) - abs(kt)]; %#ok<AGROW>
            end
        end
        RMSE(iw,ip) = sqrt(mean(err.^2, 'omitnan'));
    end
end

[minPerCol, rowIdx] = min(RMSE, [], 1, 'omitnan');
[~, colIdx] = min(minPerCol);
wBest = windowSizes(rowIdx(colIdx));
oBest = polyOrders(colIdx);

figure('Color','w');
imagesc(polyOrders, windowSizes, RMSE);
set(gca, 'YDir', 'normal');
xlabel('Polynomial order'); ylabel('Window size');
title('Curvature RMSE');
colorbar;
xticks(polyOrders); yticks(windowSizes);
hold on;
plot(oBest, wBest, 'r*', 'MarkerSize', 14, 'LineWidth', 2);

fprintf('Best: window = %d, order = %d (RMSE = %.4g)\n', wBest, oBest, RMSE(rowIdx(colIdx), colIdx));
end