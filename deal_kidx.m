function [k, idx] = deal_kidx(traj, step)
% deal_kidx  Adapter so curvature_estimation_traj2d fits the
% curvatureFcn contract of [k, idx] = curvatureFcn(x,y,z).
% traj is [x,y] (Nx2), step is the point-spacing for curvature_estimation_traj2d.
    [~, k, ~, idx] = curvature_estimation_OsculatingCircle(traj, step);
end