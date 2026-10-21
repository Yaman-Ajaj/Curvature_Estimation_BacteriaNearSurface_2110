clear; clc;

%% Settings
kThresh = 0.05;            % ignore |k| below this (noise around zero causes false flips)
minPts  = 5;               % a sign must persist this many points to count
singleSurfaceOnly = false; % true = only tracks that never visit the other surface

[fn, fp] = uigetfile('*.mat', 'Select curvature results file');
if isequal(fn,0), error('No file selected.'); end
R = load(fullfile(fp,fn));
allTrackResults = R.allTrackResults;
SF = load(R.filteredFile).AllBugs_filt;

%% Search
hits = [];   % columns: track, surface (1 = z>0, 2 = z<0), visited other surface?
for t = 1:numel(allTrackResults)
    if isempty(allTrackResults{t}) || isempty(allTrackResults{t}.sections), continue; end
    MF = SF{t};  MF = MF(~any(isnan(MF(:,2:4)),2), :);
    secs = allTrackResults{t}.sections;

    kSurf = {[], []};
    for s = 1:numel(secs)
        inSec = MF(:,1) >= secs(s).startFrame & MF(:,1) <= secs(s).endFrame;
        g = 1 + (median(MF(inSec,4)) < 0);
        kSurf{g} = [kSurf{g}; secs(s).k(:)];     % sections are in frame order
    end

    bothSurfaces = ~isempty(kSurf{1}) && ~isempty(kSurf{2});
    if singleSurfaceOnly && bothSurfaces, continue; end

    for g = 1:2
        if switchesSign(kSurf{g}, kThresh, minPts)
            hits = [hits; t, g, bothSurfaces]; %#ok<AGROW>
        end
    end
end

%% Report
fprintf('%d track/surface combinations switch curvature sign on one surface:\n', size(hits,1));
surfName = {'z > 0', 'z < 0'};
for i = 1:size(hits,1)
    fprintf('Track %5d | %s | also on other surface: %d\n', ...
        hits(i,1), surfName{hits(i,2)}, hits(i,3));
end
trackIdsToPlot = unique(hits(:,1))';   % feed this into your plotting scripts
writematrix(trackIdsToPlot(:), 'signSwitchTracks.txt');

%% Local function
function tf = switchesSign(k, kThresh, minPts)
    k = k(isfinite(k) & abs(k) > kThresh);
    if numel(k) < 2*minPts, tf = false; return; end
    s = sign(k);
    runStart = [1; find(diff(s) ~= 0) + 1];
    runLen   = diff([runStart; numel(s) + 1]);
    runSign  = s(runStart);
    runSign  = runSign(runLen >= minPts);        % keep only sustained runs
    tf = any(diff(runSign) ~= 0);
end