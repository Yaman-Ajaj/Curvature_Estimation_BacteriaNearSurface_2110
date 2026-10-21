function [k, idx] = curvature_estimation_centralDiff(traj, order)
% Estimates curvature via the parametric curvature formula
% kappa = (x'*y'' - y'*x'') / (x'^2 + y'^2)^1.5, using HIGH-ORDER
% CENTRAL DIFFERENCE formulas to estimate the first and second
% derivatives (x', y', x'', y'') at each interior point, rather than
% fitting a circle or smoothing the path.
%
% Unlike smoothTrajectorySpline-based methods, this does NOT smooth the
% data — it exactly fits a polynomial through the points in the window
% and differentiates it analytically. Higher order (larger window)
% gives a more accurate derivative estimate on clean data, but tends to
% AMPLIFY noise more than lower-order methods, since the finite-
% difference coefficients grow larger with order. Expect this method to
% do well on clean synthetic trajectories and potentially worse on
% jittered ones — that comparison is itself a useful, presentable
% result given your existing ground-truth evaluation framework.
%
% Curvature here is SIGNED (matches curvature_estimation_tangentAngle's
% convention: positive = counter-clockwise), since the parametric
% formula naturally preserves rotational sense. Take abs(k) for
% magnitude-only comparison against curvature_estimation_traj2d.
%
% Input:
%   traj  - [x, y] trajectory, n x 2. Assumes uniform frame spacing
%           (h = 1 frame) for the finite-difference formulas; if your
%           frame spacing is non-uniform this is only approximate.
%   order - which central-difference stencil to use:
%             2 - standard 3-point stencil (lowest order, least noise
%                 amplification, smallest window: step=1)
%             4 - 5-point stencil (higher accuracy, more noise
%                 amplification, window: step=2)
%           Default 4.
%
% Output:
%   k   - signed curvature at each interior point
%   idx - row indices into traj corresponding to each k value

if nargin < 2
 order = 4;
end

 n = size(traj,1);
 x = traj(:,1);
 y = traj(:,2);

switch order
case 2
 halfWidth = 1;
case 4
 halfWidth = 2;
otherwise
 error('curvature_estimation_centralDiff:unsupportedOrder', ...
 'Supported orders are 2 (3-point stencil) or 4 (5-point stencil). Got %d.', order);
end

 nPts = n - 2*halfWidth;

if nPts <= 0
 k = []; idx = [];
return;
end

 k = nan(nPts,1);
 idx = nan(nPts,1);

for cnt = 1:nPts
 i = cnt + halfWidth;
 idx(cnt) = i;

switch order
case 2
% Standard 3-point central difference (h = 1 frame spacing)
 xp = (x(i+1) - x(i-1)) / 2;
 yp = (y(i+1) - y(i-1)) / 2;
 xpp = (x(i+1) - 2*x(i) + x(i-1));
 ypp = (y(i+1) - 2*y(i) + y(i-1));

case 4
% 5-point (4th-order accurate) central difference
 xp = (-x(i+2) + 8*x(i+1) - 8*x(i-1) + x(i-2)) / 12;
 yp = (-y(i+2) + 8*y(i+1) - 8*y(i-1) + y(i-2)) / 12;
 xpp = (-x(i+2) + 16*x(i+1) - 30*x(i) + 16*x(i-1) - x(i-2)) / 12;
 ypp = (-y(i+2) + 16*y(i+1) - 30*y(i) + 16*y(i-1) - y(i-2)) / 12;
end

 speedSq = xp^2 + yp^2;

if speedSq < 1e-18
continue; % near-zero speed: heading/curvature undefined, leave as NaN
end

 k(cnt) = (xp*ypp - yp*xpp) / speedSq^1.5;
end
end