function [k, idx] = curvature_estimation_splineSmoothed_centralDiff(traj, splineP, order)

if nargin < 2
 splineP = 0.5;
end
if nargin < 3
 order = 5;
end

 trajSmoothed = smoothTrajectorySpline(traj, splineP);
 [~, k, ~, idx] = curvature_estimation_centralDiff(trajSmoothed, order);
end