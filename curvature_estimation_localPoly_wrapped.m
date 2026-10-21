function [k, idx] = curvature_estimation_localPoly_wrapped(traj, windowSize, polyOrder, t)
if nargin < 2, windowSize = 11; end
if nargin < 3, polyOrder = 3; end
if nargin < 4, t = []; end
[k, idx] = curvature_estimation_localPoly(traj, windowSize, polyOrder, t);
end