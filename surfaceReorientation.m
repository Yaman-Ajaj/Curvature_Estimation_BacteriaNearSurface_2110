% surfaceReorientation.m
% Reads a curvatureResults_*.mat file and, for every near-surface section, computes
% (all in-plane, x-y, from the filtered positions): residence time T, mean speed,
% circling rate omega two ways (mean curvature x speed, and slope of the unwrapped
% heading), swept angle, and the entry/exit headings (net displacement over nEE frames
% just outside the section) with the extra "kick" not explained by the circling.
% Sections are flagged complete (entry and exit both observed) or censored.
% The effective rotational diffusion D is estimated from the mean squared angle
% increments of the detrended heading, for surface sections and, for comparison,
% for off-surface (bulk) segments of the same tracks. Results are plotted, printed
% and saved to surfaceReorientation_<tag>.xlsx (tag from the results file name).

clear; clc;

%% Settings
umPerUnit    = 1;     % micrometers per position unit
fps          = 15;    % frames per second
nEE          = 5;     % frames outside a section used for entry/exit heading
minSecFrames = 15;    % minimum length (frames) of a segment used for D
maxLag       = 10;    % largest lag (frames) in the angle MSD fit

[fn, fp] = uigetfile('*.mat', 'Select curvature results file');
if isequal(fn,0), error('No file selected.'); end
R = load(fullfile(fp,fn));
allTrackResults = R.allTrackResults;
SF = load(R.filteredFile).AllBugs_filt;
zMin = R.surfaceZmin;  zMax = R.surfaceZmax;  tol = R.tolerance;

lags  = 1:maxLag;
sumSq = zeros(2, maxLag);  cnt = zeros(2, maxLag);   % row 1 = surface, row 2 = bulk
rows  = [];

%% Main loop
for t = 1:numel(allTrackResults)
    MF = SF{t};
    if isempty(MF), continue; end
    MF = MF(~any(isnan(MF(:,2:4)),2), :);
    if size(MF,1) < 3, continue; end
    frame = MF(:,1);  xy = MF(:,2:3);  z = MF(:,4);
    near  = (abs(z - zMax) <= tol) | (abs(z + zMin) <= tol);

    % ---- bulk segments: off-surface, frame-contiguous ----
    bKey = double(~near) .* cumsum([1; diff(frame) > 1]);
    st = find([true; diff(bKey) ~= 0]);
    en = [st(2:end) - 1; numel(bKey)];
    for r = find(bKey(st) > 0)'
        idx = st(r):en(r);
        if numel(idx) < minSecFrames, continue; end
        [tm, th] = headingSeries(frame(idx), xy(idx,:));
        if numel(tm) < 3, continue; end
        res = th - polyval(polyfit(tm, th, 1), tm);
        [sumSq(2,:), cnt(2,:)] = accumulate(res, lags, sumSq(2,:), cnt(2,:));
    end

    % ---- surface sections ----
    if isempty(allTrackResults{t}) || isempty(allTrackResults{t}.sections), continue; end
    secs = allTrackResults{t}.sections;
    for s = 1:numel(secs)
        sc = secs(s);
        inSec = frame >= sc.startFrame & frame <= sc.endFrame;
        if sum(inSec) < 4, continue; end
        [tm, th, stp] = headingSeries(frame(inSec), xy(inSec,:));
        if numel(tm) < 3, continue; end

        nFr  = sc.endFrame - sc.startFrame + 1;
        T    = nFr / fps;                                  % residence time (s)
        v    = mean(stp) * umPerUnit * fps;                % mean speed (um/s)
        surf = 1 + (median(z(inSec)) < 0);                 % 1 = z>0, 2 = z<0

        p = polyfit(tm, th, 1);
        omegaAng = p(1);                                   % rad/s, slope of unwrapped heading
        dThAng   = polyval(p, tm(end)) - polyval(p, tm(1));% fitted heading change inside section
        res = th - polyval(p, tm);
        if nFr >= minSecFrames
            [sumSq(1,:), cnt(1,:)] = accumulate(res, lags, sumSq(1,:), cnt(1,:));
        end

        kfin   = sc.k(isfinite(sc.k));
        omegaK = mean(kfin) / umPerUnit * v;               % rad/s, from curvature x speed
        dThK   = sum(kfin) * mean(stp);                    % sum(k*ds), dimensionless

        % ---- entry / exit headings from nEE frames outside the section ----
        hin = NaN;  hout = NaN;
        rIn  = ismember(frame, sc.startFrame - nEE : sc.startFrame - 1);
        rOut = ismember(frame, sc.endFrame + 1 : sc.endFrame + nEE);
        if sum(rIn) == nEE && ~any(near(rIn))
            a = find(rIn);  dd = xy(a(end),:) - xy(a(1),:);  hin = atan2(dd(2), dd(1));
        end
        if sum(rOut) == nEE && ~any(near(rOut))
            a = find(rOut); dd = xy(a(end),:) - xy(a(1),:);  hout = atan2(dd(2), dd(1));
        end
        complete = ~isnan(hin) && ~isnan(hout);
        kick = NaN;
        if complete   % exit minus entry heading, minus the circling, wrapped to (-pi, pi]
            kick = atan2(sin(hout - hin - dThAng), cos(hout - hin - dThAng));
        end

        rows = [rows; t, s, surf, complete, T, v, nFr, sc.stepReduced, ...
                omegaK, omegaAng, dThAng, dThK, hin, hout, kick]; %#ok<AGROW>
    end
end

Tb = array2table(rows, 'VariableNames', {'track','section','surface','complete', ...
    'T_s','speed_umps','nFrames','stepReduced','omega_k','omega_angle', ...
    'dTheta_angle','dTheta_k','headingIn','headingOut','kick'});

%% Rotational diffusion from the angle MSD (slope = 2D, intercept = heading noise)
msad = sumSq ./ cnt;
D = nan(1,2);  c0 = nan(1,2);
for g = 1:2
    ok = cnt(g,:) > 0;
    p = polyfit(lags(ok)/fps, msad(g,ok), 1);
    D(g) = p(1)/2;  c0(g) = p(2);
end

%% Report
isC = Tb.complete == 1;
fprintf('%d sections: %d complete, %d censored\n', height(Tb), sum(isC), sum(~isC));
for g = 1:2
    m = Tb.surface == g;
    fprintf('Surface %d: N = %d | median T = %.2f s | median omega_angle = %+.3f rad/s | median |dTheta| = %.2f rad\n', ...
        g, sum(m), median(Tb.T_s(m)), median(Tb.omega_angle(m)), median(abs(Tb.dTheta_angle(m))));
end
fprintf('Fraction of sections sweeping more than 2*pi: %.1f %%\n', 100*mean(abs(Tb.dTheta_angle) > 2*pi));
fprintf('Complete sections: kick std = %.1f deg, median |kick| = %.1f deg\n', ...
    rad2deg(std(Tb.kick(isC))), rad2deg(median(abs(Tb.kick(isC)))));
fprintf('D surface = %.4f rad^2/s | D bulk = %.4f rad^2/s | ratio = %.2f\n', D(1), D(2), D(1)/D(2));
fprintf('Heading noise floor (surface, bulk): %.3f, %.3f rad\n', sqrt(max(c0(1),0)/2), sqrt(max(c0(2),0)/2));

%% Plots
figure('Color','w'); hold on;                                   % swept angle vs residence time
scatter(Tb.T_s(isC),  Tb.dTheta_angle(isC),  18, [0.1 0.4 0.7], 'filled', 'MarkerFaceAlpha', 0.6);
scatter(Tb.T_s(~isC), Tb.dTheta_angle(~isC), 18, [0.5 0.5 0.5], 'MarkerEdgeAlpha', 0.7);
yline( 2*pi, 'k--'); yline(-2*pi, 'k--'); yline(0, 'k:');
xlabel('residence time T (s)'); ylabel('heading change inside section (rad)');
legend({'complete','censored'}); grid on; title('Swept angle vs residence time');

figure('Color','w'); hold on;                                   % omega by surface (dot plot)
for g = 1:2
    m = Tb.surface == g;
    scatter(g + (rand(sum(m),1)-0.5)*0.3, Tb.omega_angle(m), 14, 'filled', 'MarkerFaceAlpha', 0.4);
    plot(g + [-0.25 0.25], median(Tb.omega_angle(m))*[1 1], 'k-', 'LineWidth', 2);
end
yline(0, 'k:'); xlim([0.5 2.5]); set(gca, 'XTick', 1:2, 'XTickLabel', {'z > 0','z < 0'});
ylabel('\omega from heading slope (rad/s)'); grid on; title('Circling rate by surface (line = median)');

figure('Color','w'); hold on;                                   % cross-check of the two omegas
scatter(Tb.omega_angle, Tb.omega_k, 14, 'filled', 'MarkerFaceAlpha', 0.4);
lim = max(abs([Tb.omega_angle; Tb.omega_k])); plot([-lim lim], [-lim lim], 'k--');
xlabel('\omega from heading slope (rad/s)'); ylabel('\omega = mean k \cdot speed (rad/s)');
grid on; title('Cross-check of the two \omega estimates');

figure('Color','w'); hold on;                                   % angle MSD and fits
cc = {[0.85 0.33 0.10], [0.00 0.45 0.74]};  nm = {'surface','bulk'};
for g = 1:2
    plot(lags/fps, msad(g,:), 'o', 'Color', cc{g}, 'DisplayName', nm{g});
    plot(lags/fps, 2*D(g)*lags/fps + c0(g), '-', 'Color', cc{g}, 'HandleVisibility', 'off');
end
xlabel('lag (s)'); ylabel('mean squared heading increment (rad^2)'); legend; grid on;
title('Angle MSD of detrended heading (slope = 2D)');

%% Save
[~, base] = fileparts(fn);
tag = regexprep(base, '^curvatureResults_', '');
tag = regexprep(tag, '_\d{8}_\d{6}$', '');
stepTok = regexp(tag, '_step\d+', 'match', 'once');
if ~isempty(stepTok)
    tag = strrep(tag, stepTok, '');
    if contains(tag, '_sig'), tag = regexprep(tag, '_sig', [stepTok '_sig'], 'once');
    else, tag = [tag stepTok]; end
end
writetable(Tb, ['surfaceReorientation_' tag '.xlsx']);

%% Local functions
function [tm, th, stp] = headingSeries(fr, P)
% Longest frame-contiguous chunk -> unwrapped heading of one-frame steps.
% tm is the step midpoint in frames relative to the chunk start (divided by fps outside
% via the global scaling below); stp is the step length in position units.
    brk = [0; find(diff(fr) ~= 1); numel(fr)];
    [~, iL] = max(diff(brk));
    sel = (brk(iL)+1):brk(iL+1);
    fr = fr(sel);  P = P(sel,:);
    if numel(fr) < 3, tm = []; th = []; stp = []; return; end
    d   = diff(P);
    th  = unwrap(atan2(d(:,2), d(:,1)));
    stp = sqrt(sum(d.^2, 2));
    tm  = ((fr(1:end-1) + fr(2:end))/2 - fr(1)) / evalin('base', 'fps');
end

function [sumSq, cnt] = accumulate(res, lags, sumSq, cnt)
    for L = lags
        if numel(res) > L
            inc = res(1+L:end) - res(1:end-L);
            sumSq(L) = sumSq(L) + sum(inc.^2);
            cnt(L)   = cnt(L) + numel(inc);
        end
    end
end