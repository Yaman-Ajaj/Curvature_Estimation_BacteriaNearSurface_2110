function trajSmoothed = smoothTrajectorySpline(traj, sigma, t, tolMult)
% traj : N x 2 [x y]
% t    : N x 1 noise-free parameter (e.g. frame numbers). If omitted,
%        falls back to cumulative chord length.
% Fits a smoothing spline through the trajectory, parameterized by arc
% length, and evaluates it back at the SAME arc-length positions as the
% original points (no resampling) — shared helper used by any curvature
% method that wants to smooth the input geometry before differencing.
%
% Requires the Curve Fitting Toolbox (csaps).
%
% Input:
%   traj    - [x, y] trajectory, n x 2
%   splineP - csaps smoothing parameter, in [0, 1]. 0 = max smoothing,
%             1 = interpolating spline (no smoothing).
%
% Output:
%   trajSmoothed - smoothed [x, y], same size as traj, aligned 1:1 to
%                  the original points (row i of trajSmoothed
%                  corresponds to row i of traj).

% traj    : N x 2 [x y]
% sigma   : known noise std (per coordinate)
% t       : N x 1 noise-free parameter (e.g. frame numbers). If omitted
%           or empty, falls back to cumulative chord length.
% tolMult : multiplier on the theoretical tolerance N*sigma^2 (default 1)

n = size(traj,1);
if n < 3
    trajSmoothed = traj; return;
end

if nargin < 3 || isempty(t)
    df = diff(traj, 1, 1);
    t = [0; cumsum(sqrt(sum(df.^2, 2)))];
end
if nargin < 4 || isempty(tolMult)
    tolMult = 1;
end
t = t(:);

% spaps needs strictly increasing sites
[t, iu] = unique(t);
traj = traj(iu,:);
n = numel(t);
if n < 3 || t(end) == t(1)
    trajSmoothed = traj; return;
end

tol = tolMult * 2 * n * sigma^2;
sp = spaps(t, traj', tol);
trajSmoothed = fnval(sp, t)';
end