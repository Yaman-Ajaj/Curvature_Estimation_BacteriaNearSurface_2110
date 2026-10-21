function [k, idx] = curvature_estimation_localPoly(traj, windowSize, polyOrder, t)
% Local polynomial fit curvature estimator (signed).
%
% traj       : N x 2 [x y]
% windowSize : number of points in the local window (odd, >= 5)
% polyOrder  : polynomial degree (must be < windowSize)
% t          : N x 1 parameter (e.g. frame numbers). If omitted, uses
%              point index 1:N (uniform spacing assumed).

if nargin < 2, windowSize = 11; end
if nargin < 3, polyOrder = 3; end
if nargin < 4 || isempty(t)
    t = (1:size(traj,1))';
end
t = t(:);

if mod(windowSize,2) == 0
    windowSize = windowSize + 1;
    warning('curvature_estimation_localPoly:windowMustBeOdd', ...
        'windowSize must be odd; using %d.', windowSize);
end
if polyOrder >= windowSize
    error('curvature_estimation_localPoly:orderTooHigh', ...
        'polyOrder (%d) must be less than windowSize (%d).', polyOrder, windowSize);
end

n = size(traj,1);
halfWindow = floor(windowSize/2);
nPts = n - 2*halfWindow;
if nPts <= 0
    k = []; idx = [];
    return;
end

x = traj(:,1);
y = traj(:,2);

k = nan(nPts,1);
idx = nan(nPts,1);

for cnt = 1:nPts
    i = cnt + halfWindow;
    ind = (i-halfWindow):(i+halfWindow);
    tau = t(ind) - t(i);   % local coordinate centered at point i

    px = polyfit(tau, x(ind), polyOrder);
    py = polyfit(tau, y(ind), polyOrder);

    dpx = polyder(px);   dpy = polyder(py);
    ddpx = polyder(dpx); ddpy = polyder(dpy);

    dx = polyval(dpx, 0);
    dy = polyval(dpy, 0);
    ddx = polyval(ddpx, 0);
    ddy = polyval(ddpy, 0);

    denom = (dx^2 + dy^2)^1.5;
    if denom > eps
        k(cnt) = (dx*ddy - dy*ddx) / denom;   % signed, no abs()
    end
    idx(cnt) = i;
end
end