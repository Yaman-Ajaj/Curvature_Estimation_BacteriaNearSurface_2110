function [k, idx, cc] = curvature_estimation_splineDerivative(traj, sigma, t, tolMult, m)
% Computes signed curvature analytically from spline derivatives.
%
% traj    : N x 2 [x y]
% sigma   : known noise std (per coordinate)
% t       : N x 1 noise-free parameter (e.g. frame numbers). If omitted
%           or empty, falls back to cumulative chord length.
% tolMult : multiplier on the theoretical tolerance (default 1)
%
% OUTPUT:
%   k   : signed curvature at each returned point
%   idx : row indices into traj (all points; no edge trimming needed)
%   cc  : [NaN NaN] placeholder for compatibility with circle-fit outputs
%         (no explicit center is computed with this method)

if nargin < 3, t = []; end
if nargin < 4, tolMult = 1; end
if nargin < 5, m = 3; end

n = size(traj,1);
if n < 4
    k = []; idx = []; cc = [];
    return;
end

if isempty(t)
    df = diff(traj, 1, 1);
    t = [0; cumsum(sqrt(sum(df.^2, 2)))];
end
t = t(:);

[t, iu] = unique(t);
traj = traj(iu,:);
n = numel(t);
if n < 4 || t(end) == t(1)
    k = []; idx = []; cc = [];
    return;
end

tol = tolMult * 2 * n * sigma^2;
sp = spaps(t, traj', tol, [],m);

d1 = fnder(sp, 1);   % velocity
d2 = fnder(sp, 2);   % acceleration

v  = fnval(d1, t)';  % N x 2: [x' y']
a  = fnval(d2, t)';  % N x 2: [x'' y'']

num   = v(:,1).*a(:,2) - v(:,2).*a(:,1);
denom = (v(:,1).^2 + v(:,2).^2).^1.5;

k = nan(n,1);
valid = denom > eps;
k(valid) = num(valid) ./ denom(valid);

idx = iu; % maps back to positions in the (unique-filtered) input traj
cc  = nan(n,2);
end