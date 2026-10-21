function [k, idx] = curvatureEstimationTangentAngle(traj, step)
% Estimates curvature at each interior point via the derivative of the
% tangent (heading) angle with respect to arc length: k = d(alpha)/ds.
%
% Uses a step-offset window (points i-step and i+step relative to i),
% matching the same baseline convention as curvature_estimation_traj2d,
% so the two methods can be compared at the same spatial baseline
% rather than one using immediate neighbors and the other a wider
% window. Widening step reduces noise sensitivity but also blurs any
% real curvature variation happening within that span — for that
% reason, prefer keeping step small and applying smoothCurvature to
% the output afterward if less noisy results are needed, rather than
% relying on step alone to suppress noise.
%
% Input:
%   traj - [x, y] trajectory, n x 2
%   step - index offset (points before/after) used to build the two
%          step vectors at each point. Default 1 (immediate neighbors,
%          matching the original un-parameterized version).
%
% Output:
%   k   - signed curvature at each interior point (n-2*step x 1).
%         Sign reflects rotational sense (positive = counter-clockwise,
%         negative = clockwise), unlike curvature_estimation_traj2d's
%         k, which is always positive. Take abs(k) for magnitude-only
%         comparison against that method.
%         Points where either step has ~zero displacement (heading
%         undefined) are reported as NaN, not skipped — idx stays
%         gap-free and aligned with k.
%   idx - row indices into traj corresponding to each k value.

if nargin < 2
 step = 1;
end

 n = size(traj,1);
 nPts = n - 2*step;

if nPts <= 0
 k = []; idx = [];
return;
end

 k = nan(nPts,1);
 idx = nan(nPts,1);

 minDisplacement = 1e-9;

for cnt = 1:nPts
 i = cnt + step; % interior point index into traj

 v1 = traj(i,:) - traj(i-step,:); % incoming step, spanning `step` points
 v2 = traj(i+step,:) - traj(i,:); % outgoing step, spanning `step` points

 ds1 = norm(v1);
 ds2 = norm(v2);

 idx(cnt) = i;

if ds1 < minDisplacement || ds2 < minDisplacement
continue;
end

 alpha1 = atan2(v1(2), v1(1));
 alpha2 = atan2(v2(2), v2(1));

 dAlpha = atan2(sin(alpha2 - alpha1), cos(alpha2 - alpha1));

 ds = (ds1 + ds2) / 2;

 k(cnt) = dAlpha / ds;
end
end