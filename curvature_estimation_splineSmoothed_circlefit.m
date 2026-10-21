function [k, idx] = curvature_estimation_splineSmoothed_circlefit(traj, sigma, step, t, tolMult)
% Same input-smoothing idea as curvature_estimation_splineSmoothed, but
% runs the smoothed path through the osculating-circle method
% (curvature_estimation_traj2d) instead of tangent-angle differencing —
% isolating "does input-smoothing help" as a question independent of
% which curvature formula is used afterward.
%
% Returns [k, idx] only (r and cc are discarded), matching the
% curvatureFcn contract used directly by the main loop for every
% method except the raw circle-fit method, which uses the separate
% deal_kidx adapter instead.

if nargin < 2, sigma = 0.2; end
if nargin < 3, step = 1; end
if nargin < 4, t = []; end
if nargin < 5, tolMult = 1; end

trajSmoothed = smoothTrajectorySpline(traj, sigma, t, tolMult);
[~, k, ~, idx] = curvature_estimation_OsculatingCircle(trajSmoothed, step);
end