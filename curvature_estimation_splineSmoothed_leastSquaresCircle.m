function [k, idx] = curvature_estimation_splineSmoothed_leastSquaresCircle(traj, sigma, windowSize, t, tolMult)
% Same input-smoothing idea as curvature_estimation_splineSmoothed and
% curvature_estimation_splineSmoothed_circlefit, but runs the smoothed
% path through the least-squares circle fit
% (curvature_estimation_leastSquaresCircle) instead of tangent-angle
% differencing or an exact 3-point circle fit — isolating "does
% input-smoothing help" as a question independent of which curvature
% formula is used afterward, now across all three formulas.
%
% Returns [k, idx] only (r and cc are discarded), matching the
% curvatureFcn contract, same pattern as
% curvature_estimation_splineSmoothed_circlefit.

    
if nargin < 2, sigma = 0.2; end
if nargin < 3, windowSize = 5; end
if nargin < 4, t = []; end
if nargin < 5, tolMult = 1; end

trajSmoothed = smoothTrajectorySpline(traj, sigma, t, tolMult);
[~, k, ~, idx] = curvature_estimation_leastSquaresCircle(trajSmoothed, windowSize);
end