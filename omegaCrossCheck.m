% omegaCrossCheck.m
% Cross-check: fit log|omega| vs log(speed), omega = turning rate of the heading
% angle (from filtered x,y positions, no curvature estimator).
% Expectation: slope_omega = 1 + b, where b is the log k vs log v exponent.
% Requires the Statistics Toolbox (robustfit).

clear; rng(1);

%% Settings (keep identical to curvatureVsSpeed.m)
framesPerSec = 15;
umPerUnit    = 1;
useAbs       = false;  % true: mean|dtheta| and mean|k|; false: |mean dtheta| and |mean k|
minPts       = 5;      % min finite k values per section
cleanOnly    = false;
minFrames    = 15;
nBoot        = 1000;

[fn, fp] = uigetfile('*.mat', 'Select curvature results file');
if isequal(fn,0), error('No file selected.'); end
R = load(fullfile(fp,fn));
allTrackResults = R.allTrackResults;
SF = load(R.filteredFile).AllBugs_filt;

%% Per-section: [track, |k|, |omega| (rad/s), speed]
rows = [];
for t = 1:numel(allTrackResults)
    if isempty(allTrackResults{t}) || isempty(allTrackResults{t}.sections), continue; end
    MF = SF{t};  MF = MF(~any(isnan(MF(:,2:4)),2), :);
    secs = allTrackResults{t}.sections;

    for s = 1:numel(secs)
        k = secs(s).k(:);  k = k(isfinite(k)) / umPerUnit;
        if numel(k) < minPts, continue; end

        inSec = MF(:,1) >= secs(s).startFrame & MF(:,1) <= secs(s).endFrame;
        nIn = sum(inSec);
        if nIn < 3, continue; end
        if cleanOnly && (secs(s).stepReduced || nIn < minFrames), continue; end

        fr  = MF(inSec,1);
        stp = diff(MF(inSec,2:3));                 % in-plane step vectors
        ok  = diff(fr) == 1;                       % never use steps across gaps
        if sum(ok) < 2, continue; end

        % speed (same definition as curvatureVsSpeed.m)
        spd = mean(sqrt(sum(stp(ok,:).^2, 2))) * umPerUnit * framesPerSec;

        % heading change between consecutive gap-free steps (wrapped to [-pi, pi])
        th    = atan2(stp(:,2), stp(:,1));
        dth   = angle(exp(1i*diff(th)));
        dth   = dth(ok(1:end-1) & ok(2:end));
        if isempty(dth), continue; end

        if useAbs
            om = mean(abs(dth));  kM = mean(abs(k));
        else
            om = abs(mean(dth));  kM = abs(mean(k));
        end
        rows = [rows; t, kM, om*framesPerSec, spd]; %#ok<AGROW>
    end
end

%% Fits
good = rows(:,2) > 0 & rows(:,3) > 0 & rows(:,4) > 0;
rows = rows(good,:);
trk = rows(:,1);  k = rows(:,2);  om = rows(:,3);  v = rows(:,4);
lv = log(v);

pK  = robustfit(lv, log(k));   ciK = bootRobustB(lv, log(k),  trk, nBoot);
pO  = robustfit(lv, log(om));  ciO = bootRobustB(lv, log(om), trk, nBoot);
olsO = polyfit(lv, log(om), 1);
rho = corr(v, om, 'Type', 'Spearman');

fprintf('\nN sections = %d\n', numel(v));
fprintf('b (log k vs log v, robust)           = %+.2f  (95%% CI %+.2f to %+.2f)\n', pK(2), ciK);
fprintf('1 + b (predicted omega slope)        = %+.2f  (95%% CI %+.2f to %+.2f)\n', 1+pK(2), 1+ciK);
fprintf('omega slope, direct (robust)         = %+.2f  (95%% CI %+.2f to %+.2f)\n', pO(2), ciO);
fprintf('omega slope, direct (OLS)            = %+.2f\n', olsO(1));
fprintf('Spearman rho (omega vs speed)        = %+.3f\n', rho);
fprintf('Difference (direct - predicted)      = %+.2f\n', pO(2) - (1+pK(2)));

%% Plot
vv = linspace(min(v), max(v), 100);
figure('Color','w');
loglog(v, om, '.', 'Color', [0.7 0.7 0.7]); hold on;
loglog(vv, exp(pO(1)) * vv.^pO(2), 'r-', 'LineWidth', 1.5);
loglog(vv, exp(pO(1)) * vv.^(1+pK(2)), 'b--', 'LineWidth', 1.2);  % predicted slope, same intercept
xlabel('mean speed of section'); ylabel('|\omega| of section (rad/s)'); grid on;
legend({'sections','direct robust fit','slope predicted by 1 + b'}, 'Location', 'southwest');
title('Angular velocity vs speed (log-log)');

%% Local function
function ci = bootRobustB(lv, ly, trk, nBoot)
% By-track bootstrap 95% CI of the robust slope.
[~, ~, g] = unique(trk);
members = accumarray(g, (1:numel(g))', [], @(x){x});
nT = numel(members);
bB = nan(nBoot,1);
for i = 1:nBoot
    j = vertcat(members{randi(nT, nT, 1)});
    pp = robustfit(lv(j), ly(j));
    bB(i) = pp(2);
end
ci = prctile(bB, [2.5 97.5]);
end