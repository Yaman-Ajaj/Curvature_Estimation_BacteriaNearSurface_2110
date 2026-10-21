function [k, idx] = curvature_estimation_splineSmoothed(traj, sigma, step, t, tolMult)
if nargin < 2, sigma = 0.2; end
if nargin < 3, step = 1; end
if nargin < 4, t = []; end
if nargin < 5, tolMult = 1; end

trajSmoothed = smoothTrajectorySpline(traj, sigma, t, tolMult);
[k, idx] = curvatureEstimationTangentAngle(trajSmoothed, step);
end