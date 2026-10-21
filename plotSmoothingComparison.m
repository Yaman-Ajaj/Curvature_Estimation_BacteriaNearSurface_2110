clear; clc;

%% Settings
dataFile = 'Ref0722_unfiltered.mat';   % <-- file containing S.Speeds
%trackIds = [4383 3461 1981];           % tracks to plot
rng(1)
nTracksAvailable = numel(load(dataFile).AllBugs_unfilt);
trackIds = randperm(nTracksAvailable, 7);

sigma    = 0.25;
tolMult  = 1.29;                       % use the value you analyse with
surfaceZmin = 35.9535; surfaceZmax = 35.8652; tolerance = 3.5;
minLen = 7;                            % skip shorter runs

S = load(dataFile).AllBugs_unfilt;
SF = load('Ref0722_filtered.mat').AllBugs_filt;

 for trackIdx = trackIds
        M  = S{trackIdx};   MF = SF{trackIdx};
    if isempty(M) || isempty(MF), continue; end
    M  = M(~any(isnan(M(:,2:4)),2), :);
    MF = MF(~any(isnan(MF(:,2:4)),2), :);
    if isempty(M) || isempty(MF), continue; end
    frame  = M(:,1);  x = M(:,2);  y = M(:,3);   % unfiltered: what gets smoothed
    frameF = MF(:,1); zF = MF(:,4);              % filtered: used only for the z test

    nearSurface = (abs(zF - surfaceZmax) <= tolerance) | (abs(zF + surfaceZmin) <= tolerance);
    runKey = double(nearSurface) .* cumsum([1; diff(frameF) > 1]);
    starts = find([true; diff(runKey) ~= 0]);
    ends   = [starts(2:end) - 1; numel(runKey)];
    keep   = runKey(starts) > 0;
    starts = starts(keep);  ends = ends(keep);

    % map each filtered run to rows of the unfiltered data, split at unfiltered frame gaps
    runsU = {};
    for r = 1:numel(starts)
        rowsU = find(ismember(frame, frameF(starts(r):ends(r))));
        if isempty(rowsU), continue; end
        sub = cumsum([1; diff(frame(rowsU)) > 1]);
        for q = unique(sub)'
            runsU{end+1} = rowsU(sub == q); %#ok<AGROW>
        end
    end

    figure('Color','w','Name',sprintf('Track %d',trackIdx));
    axXY = subplot(1,2,1); hold(axXY,'on'); axis(axXY,'equal'); grid(axXY,'on');
    axT  = subplot(1,2,2); hold(axT,'on');  grid(axT,'on');
    nPlotted = 0;

    for r = 1:numel(runsU)
        idx = runsU{r};
        if numel(idx) < minLen, continue; end
        traj = [x(idx) y(idx)];
        t = frame(idx);
        sm = smoothTrajectorySpline(traj, sigma, t, tolMult);

        plot(axXY, traj(:,1), traj(:,2), '.-', 'Color',[0.7 0.7 0.7], 'MarkerSize',8);
        plot(axXY, sm(:,1), sm(:,2), 'r.-', 'LineWidth',1.3, 'MarkerSize',8);

        plot(axT, t, traj(:,1), '.', 'Color',[0.6 0.6 1]);
        plot(axT, t, sm(:,1), 'b-', 'LineWidth',1.3);
        plot(axT, t, traj(:,2), '.', 'Color',[1 0.6 0.6]);
        plot(axT, t, sm(:,2), 'r-', 'LineWidth',1.3);
        nPlotted = nPlotted + 1;
    end

    xlabel(axXY,'x'); ylabel(axXY,'y');
    title(axXY, 'Track : raw (grey) vs smoothed (red)');
    xlabel(axT,'frame'); ylabel(axT,'position');
    title(axT, 'x (blue) and y (red) vs frame: dots raw, lines smoothed');
    sgtitle(sprintf('Track %d  |  \\sigma = %g,  tolMult = %g', trackIdx, sigma, tolMult));
    if nPlotted == 0, close(gcf); fprintf('Track %d: no run >= %d points\n', trackIdx, minLen); end
end