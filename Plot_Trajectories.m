%Plots the trajectories that have the indicies in the list, each trajectory
%is in a different plot

% Load data
data = load('Ref0722_filtered.mat');
S = data.AllBugs_filt;


kList = [8];
%kList = randperm(10509, 10);

for i = 1:length(kList)
    k = kList(i);

    M = data.AllBugs_filt{k};

    % Validate
    if isempty(M) || size(M,2) < 4
        warning(['Bacterium ', num2str(k), ' skipped (invalid data)']);
        continue;
    end

    % Extract coordinates
    x = M(:,2);
    y = M(:,3);
    z = M(:,4);
    frame = M(:,1);

    figure; % NEW figure for each bacterium

    scatter3(x, y, z, 20, frame, 'filled');
    colorbar;

    xlabel('x');
    ylabel('y');
    zlabel('z');
    title(['3D trajectory (colored by time), bacterium ', num2str(k)]);

    grid on;
    axis equal;
    view(3);
end
