% Load data
%data = load('Fake_circular_trajectories_noisy.mat');
data = load('Fake_synthetic_trajectories_Noise.mat');
trajectories = data.trajectories;

rng(2)
kList = randperm(40, 10);
%kList = [4];
for i = 1:length(kList)
    k = kList(i);
    M = trajectories{k};

    % Validate
    if isempty(M) || size(M,2) < 4
        warning(['Trajectory ', num2str(k), ' skipped (invalid data)']);
        continue;
    end

    % Extract coordinates
    x = M(:,2);
    y = M(:,3);
    z = M(:,4);
    frame = M(:,1);

    figure; % NEW figure for each trajectory
    scatter3(x, y, z, 20, frame, 'filled');
    colorbar;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(['3D trajectory (colored by time), trajectory ', num2str(k)]);
    grid on;
    axis equal;
    view(3);
end