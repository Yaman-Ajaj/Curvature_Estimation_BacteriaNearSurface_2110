function [r,k,cc,idx] = curvature_estimation_OsculatingCircle(traj, step, r_max)
% Estimates curvature at each interior trajectory point via the 3-point
% circumscribed-circle method (perpendicular bisectors of consecutive
% chords). Points step samples apart are used to reduce noise sensitivity.
%
% traj : N x 2 [x y]
% step : spacing (in samples) between the center point and its neighbors
% on each side (default 1)
% r_max : radius reported for near-collinear points, where the bisectors
% are (near-)parallel and no finite circumcenter exists (default 1e6)
%
% OUTPUT:
% r : circumradius at each interior point (r_max if near-collinear)
% k : curvature = 1/r at each interior point
% cc : [x y] circumcenter at each interior point ([NaN NaN] if near-collinear)
% idx : row indices into traj corresponding to each output (edges trimmed)

if nargin < 2
 step = 1;
end
if nargin < 3
 r_max = 1e6; % large finite radius reported for near-collinear (straight-swimming) points
end
n = size(traj,1);
nPts = n - 2*step;
if nPts <= 0
 r = []; k = []; cc = []; idx = [];
return;
end
r = nan(nPts,1);
k = nan(nPts,1);
cc = nan(nPts,2);
idx = nan(nPts,1);
cnt = 0;
for i = 1+step : n-step %Only interior points are used, neglecting the two points at the edges
 ixL = i-step;
 ixC = i;
 ixR = i+step;
 pi_kb = traj(ixL,:);
 pi_kf = traj(ixR,:);
 pi_k = traj(ixC,:);
 [a1,b1,c1] = points_bisect_line_imp_2d(pi_k,pi_kb);
 [a2,b2,c2] = points_bisect_line_imp_2d(pi_k,pi_kf);
 [ival,center] = lines_imp_int_2d(a1,b1,c1,a2,b2,c2);
 cnt = cnt + 1;
if ival == 1
    radius = norm(pi_k - center);
    r(cnt) = radius; k(cnt) = 1/radius; cc(cnt,:) = center;
elseif ival == 0
    % parallel bisectors: collinear points, straight segment
    r(cnt) = r_max; k(cnt) = 1/r_max; cc(cnt,:) = [NaN NaN];
else
    % ival = -1, -2, 2: coincident points, curvature undefined
    r(cnt) = NaN; k(cnt) = NaN; cc(cnt,:) = [NaN NaN];
end
 idx(cnt) = ixC;
end
% cnt should always equal nPts since every iteration now writes a value;
% this trim is just a safety net in case that assumption is ever violated
r = r(1:cnt);
k = k(1:cnt);
cc = cc(1:cnt,:);
idx = idx(1:cnt);
end

function [ a, b, c ] = points_bisect_line_imp_2d ( p1, p2 )
  dim_num = 2;
  a = p1(1) - p2(1);
  b = p1(2) - p2(2);
  c = - 0.5 * ( sum ( p1(1:dim_num).^2 ) - sum ( p2(1:dim_num).^2 ) );
end

function [ ival, p ] = lines_imp_int_2d ( a1, b1, c1, a2, b2, c2 )
  dim_num = 2;
if ( a1 == 0.0 && b1 == 0.0 )
    ival = -1;
    p = [];
return
elseif ( a2 == 0.0 && b2 == 0.0 )
    ival = -2;
    p = [];
return
end
  a(1,1) = a1;
  a(1,2) = b1;
  a(1,3) = -c1;
  a(2,1) = a2;
  a(2,2) = b2;
  a(2,3) = -c2;
  [ a, info ] = r8mat_solve ( 2, 1, a );
if ( info == 0 )
    ival = 1;
    p(1:dim_num) = a(1:dim_num,3);
else
    ival = 0;
    p = [];
if ( a1 == 0.0 )
if ( b2 * c1 == c2 * b1 )
        ival = 2;
        p(1:dim_num) = [ 0.0, - c1 / b1 ];
end
else
if ( a2 * c1 == c2 * a1 )
        ival = 2;
if ( abs ( a1 ) < abs ( b1 ) )
          p(1:dim_num) = [ 0.0, - c1 / b1 ];
else
          p(1:dim_num) = [ - c1 / a1, 0.0 ];
end
end
end
end
end

function [ a, info ] = r8mat_solve ( n, nrhs, a )
  info = 0;
for j = 1 : n
    ipivot = j;
    apivot = a(j,j);
for i = j+1 : n
if ( abs ( apivot ) < abs ( a(i,j) ) )
        apivot = a(i,j);
        ipivot = i;
end
end
if ( apivot == 0.0 )
      info = j;
return;
end
    temp               = a(ipivot,1:n+nrhs);
    a(ipivot,1:n+nrhs) = a(j,     1:n+nrhs);
    a(j,     1:n+nrhs) = temp;
    a(j,j) = 1.0;
    a(j,j+1:n+nrhs) = a(j,j+1:n+nrhs) / apivot;
for i = 1 : n
if ( i ~= j )
        factor = a(i,j);
        a(i,j) = 0.0;
        a(i,j+1:n+nrhs) = a(i,j+1:n+nrhs) - factor * a(j,j+1:n+nrhs);
end
end
end
end