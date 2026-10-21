function [k, idx] = curvature_estimation_splineSmoothed_derivative(traj, sigma, t, tolMult)
if nargin < 2, sigma = 0.2; end
if nargin < 3, t = []; end
if nargin < 4, tolMult = 1; end

[k, idx] = curvature_estimation_splineDerivative(traj, sigma, t, tolMult);
end