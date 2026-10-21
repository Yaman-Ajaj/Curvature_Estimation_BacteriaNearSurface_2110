function [r, k, cc, idx] = curvature_estimation_leastSquaresCircle(traj, windowSize, r_max)
% Estimates curvature at each interior point by fitting a circle via
% ALGEBRAIC LEAST SQUARES (Kasa method) through a window of
% (2*windowHalf+1) points centered on that point, rather than an exact
% fit through only 3 points. Unlike curvature_estimation_traj2d, this
% is an OVER-DETERMINED fit (more points than the 3 circle parameters),
% so noise in individual points is averaged out rather than forcing the
% circle through every point exactly.
%
% windowSize is a window-size knob in the same spirit as curvatureStep
% or smoothWindow elsewhere in this pipeline: larger windows average
% out more noise but blur real curvature variation happening within
% the window. Tune empirically via ground-truth comparison, the same
% way as your other methods.
%
% Input:
%   traj       - [x, y] trajectory, n x 2
%   windowSize - number of points used in each local circle fit.
%                Must be odd and >= 5 (an even-numbered or smaller
%                window is auto-adjusted with a warning). Default 5.
%   r_max      - placeholder radius reported when a window's points are
%                too close to collinear for a stable fit (mirrors
%                curvature_estimation_traj2d's r_max convention).
%                Default 1e6.
%
% Output:
%   r   - fitted radius at each interior point (n-2*windowHalf x 1)
%   k   - curvature (1/r) at each point, always positive (unsigned),
%         matching curvature_estimation_traj2d's convention
%   cc  - fitted circle center [x,y] at each point (n-2*windowHalf x 2)
%   idx - row indices into traj corresponding to each output, same
%         convention as curvature_estimation_traj2d

if nargin < 2
 windowSize = 5;
end
if nargin < 3
 r_max = 1e6;
end

if windowSize < 5
 warning('curvature_estimation_leastSquaresCircle:windowTooSmall', ...
 'windowSize must be >= 5 for an over-determined fit; using 5.');
 windowSize = 5;
end
if mod(windowSize, 2) == 0
 windowSize = windowSize + 1;
 warning('curvature_estimation_leastSquaresCircle:windowMustBeOdd', ...
 'windowSize must be odd (centered window); using %d.', windowSize);
end

 windowHalf = (windowSize - 1) / 2;
 n = size(traj,1);
 nPts = n - 2*windowHalf;

if nPts <= 0
 r = []; k = []; cc = []; idx = [];
return;
end

 r = nan(nPts,1);
 k = nan(nPts,1);
 cc = nan(nPts,2);
 idx = nan(nPts,1);

for cnt = 1:nPts
 i = cnt + windowHalf;
 windowIdx = (i - windowHalf) : (i + windowHalf);
 pts = traj(windowIdx, :);

 idx(cnt) = i;

 [center, radius, isValid] = fitCircleKasa(pts);

if isValid
 r(cnt) = radius;
 k(cnt) = 1/radius;
 cc(cnt,:) = center;
% else
% % Points too close to collinear for a stable fit — mirror
% % curvature_estimation_traj2d's near-straight handling
%  r(cnt) = r_max;
%  k(cnt) = 1/r_max;
%  cc(cnt,:) = [NaN NaN];
end
end
end

function [center, radius, isValid] = fitCircleKasa(pts, r_max)
% Algebraic least-squares circle fit (Kasa method).
%
% Model: x^2 + y^2 + D*x + E*y + F = 0  (general circle equation,
% expanded), solved as a linear least-squares problem for [D, E, F]
% given the points. Center = (-D/2, -E/2), radius from center + F.
%
% This is fast (one linear solve, no iteration) but has a known mild
% bias toward overestimating radius when points span only a small arc
% of the circle rather than going most of the way around — worth
% keeping in mind when interpreting results on very gentle arcs.
    if nargin < 2, r_max = 1e6; end
    mu = mean(pts, 1);
    x = pts(:,1) - mu(1);
    y = pts(:,2) - mu(2);

    A = [x, y, ones(numel(x),1)];
    if rank(A) < 3
        center = [NaN NaN]; radius = NaN; isValid = false; return;      % collinear/coincident
    end

    sol = A \ (-(x.^2 + y.^2));
    cx = -sol(1)/2;  cy = -sol(2)/2;
    radiusSq = cx^2 + cy^2 - sol(3);

    if radiusSq <= 0
        center = [NaN NaN]; radius = NaN; isValid = false; return;      % failed fit
    elseif sqrt(radiusSq) > r_max
        center = [NaN NaN]; radius = Inf; isValid = true; return;       % straight, k = 0
    end

    center = [cx cy] + mu;      % back to original coordinates
    radius = sqrt(radiusSq);
    isValid = true;
end