clear; clc;
rng(2)

%% Adjustable parameters
nPerType = 10;
framesRange = [50, 200];
radiusRange = [5, 20];
semiAxisRange = [5, 20];
arcSpanRangeDeg = [45, 150];

% SpatiallyVarying parameters (gentle, smooth curvature drift)
svTSpanRange = [15, 30];
svFreqRange  = [0.3, 1.2];
svAmpRange   = [0.05, 0.2];

% FastVarying parameters (short-wavelength curvature oscillation)
fvCyclesRange = [4, 8];      % number of full curvature oscillations per track
fvAmpRange    = [0.05, 0.15];
fvBaseRange   = [0.02, 0.08];

centerX = 0;
centerY = 0;
zValue = 36;
startFrame = 1;

outputFileNameClean = 'Fake_synthetic_trajectories_Clean.mat';
outputFileNameNoisy = 'Fake_synthetic_trajectories_Noise.mat';
sigmaTotal = 0.25;                                 % effective noise estimated from real data
quantStep  = 0.325;                                % pixel grid of the imaging
jitterStdNoisy = sqrt(sigmaTotal^2 - quantStep^2/12);   % ≈ 0.22

shapeTypes = {'Circular', 'GentleArc', 'Elliptical', 'SpatiallyVarying', %'FastVarying'
    };
nTotal = nPerType * numel(shapeTypes);

%% --- Pass 1: generate ALL shape parameters once, before any jitter ---
paramFrameCounts = zeros(nTotal, 1);
paramRadii = nan(nTotal, 1);
paramSemiAxes = nan(nTotal, 2);
paramArcSpanDeg = nan(nTotal, 1);
paramType = cell(nTotal, 1);

paramSV = struct('tSpan', nan(nTotal,1), 'freq1', nan(nTotal,1), 'freq2', nan(nTotal,1), ...
                  'amp1', nan(nTotal,1), 'amp2', nan(nTotal,1));

paramFV = struct('nCycles', nan(nTotal,1), 'kAmp', nan(nTotal,1), 'kBase', nan(nTotal,1));

idx = 0;
for i = 1:nPerType
    idx = idx + 1;
    paramType{idx} = 'Circular';
    paramFrameCounts(idx) = randi(framesRange);
    paramRadii(idx) = radiusRange(1) + (radiusRange(2) - radiusRange(1)) * rand();
    paramArcSpanDeg(idx) = 360;
end
for i = 1:nPerType
    idx = idx + 1;
    paramType{idx} = 'GentleArc';
    paramFrameCounts(idx) = randi(framesRange);
    paramRadii(idx) = radiusRange(1) + (radiusRange(2) - radiusRange(1)) * rand();
    paramArcSpanDeg(idx) = arcSpanRangeDeg(1) + (arcSpanRangeDeg(2) - arcSpanRangeDeg(1)) * rand();
end
for i = 1:nPerType
    idx = idx + 1;
    paramType{idx} = 'Elliptical';
    paramFrameCounts(idx) = randi(framesRange);
    axis1 = semiAxisRange(1) + (semiAxisRange(2) - semiAxisRange(1)) * rand();
    axis2 = semiAxisRange(1) + (semiAxisRange(2) - semiAxisRange(1)) * rand();
    paramSemiAxes(idx, :) = [max(axis1, axis2), min(axis1, axis2)];
end
for i = 1:nPerType
    idx = idx + 1;
    paramType{idx} = 'SpatiallyVarying';
    paramFrameCounts(idx) = randi(framesRange);
    paramSV.tSpan(idx) = svTSpanRange(1) + diff(svTSpanRange)*rand();
    paramSV.freq1(idx) = svFreqRange(1) + diff(svFreqRange)*rand();
    paramSV.freq2(idx) = svFreqRange(1) + diff(svFreqRange)*rand();
    paramSV.amp1(idx)  = svAmpRange(1) + diff(svAmpRange)*rand();
    paramSV.amp2(idx)  = svAmpRange(1) + diff(svAmpRange)*rand();
end
for i = 1:nPerType
    idx = idx + 1;
    paramType{idx} = 'FastVarying';
    paramFrameCounts(idx) = randi(framesRange);
    paramFV.nCycles(idx) = randi(fvCyclesRange);
    paramFV.kAmp(idx)    = fvAmpRange(1) + diff(fvAmpRange)*rand();
    paramFV.kBase(idx)   = fvBaseRange(1) + diff(fvBaseRange)*rand();
end

%% --- Pass 2: build both clean and noisy trajectory sets from the SAME parameters ---
[trajectoriesClean, trajectoryType, trajectoryRadii, trajectorySemiAxes, ...
    trajectoryArcSpanDeg, trajectoryFrameCounts, trajectoryTheta, trajectorySV, ...
    trajectoryKappaTrue, trajectoryFV] = buildTrajectorySet( ...
    paramType, paramFrameCounts, paramRadii, paramSemiAxes, paramArcSpanDeg, paramSV, paramFV, ...
    centerX, centerY, zValue, startFrame, 0,0); % jitterStd = 0

trajectories = trajectoriesClean;
save(outputFileNameClean, 'trajectories', 'trajectoryType', 'trajectoryRadii', ...
    'trajectorySemiAxes', 'trajectoryArcSpanDeg', 'trajectoryFrameCounts', 'trajectoryTheta', ...
    'trajectorySV', 'trajectoryKappaTrue', 'trajectoryFV');
disp(['Saved file: ', outputFileNameClean]);

[trajectoriesNoisy, trajectoryType, trajectoryRadii, trajectorySemiAxes, ...
    trajectoryArcSpanDeg, trajectoryFrameCounts, trajectoryTheta, trajectorySV, ...
    trajectoryKappaTrue, trajectoryFV] = buildTrajectorySet( ...
    paramType, paramFrameCounts, paramRadii, paramSemiAxes, paramArcSpanDeg, paramSV, paramFV, ...
    centerX, centerY, zValue, startFrame, jitterStdNoisy, quantStep);

trajectories = trajectoriesNoisy;
save(outputFileNameNoisy, 'trajectories', 'trajectoryType', 'trajectoryRadii', ...
    'trajectorySemiAxes', 'trajectoryArcSpanDeg', 'trajectoryFrameCounts', 'trajectoryTheta', ...
    'trajectorySV', 'trajectoryKappaTrue', 'trajectoryFV');
disp(['Saved file: ', outputFileNameNoisy]);

%% --- Helper: builds one full trajectory set from fixed shape parameters ---
function [trajectories, trajectoryType, trajectoryRadii, trajectorySemiAxes, ...
    trajectoryArcSpanDeg, trajectoryFrameCounts, trajectoryTheta, trajectorySV, ...
    trajectoryKappaTrue, trajectoryFV] = buildTrajectorySet( ...
    paramType, paramFrameCounts, paramRadii, paramSemiAxes, paramArcSpanDeg, paramSV, paramFV, ...
    centerX, centerY, zValue, startFrame, jitterStd, quantStep)

    nTotal = numel(paramType);
    trajectories = cell(nTotal, 1);
    trajectoryTheta = cell(nTotal, 1);
    trajectoryKappaTrue = cell(nTotal, 1);
    trajectorySV = paramSV;
    trajectoryFV = paramFV;
    trajectoryType = paramType;
    trajectoryRadii = paramRadii;
    trajectorySemiAxes = paramSemiAxes;
    trajectoryArcSpanDeg = paramArcSpanDeg;
    trajectoryFrameCounts = paramFrameCounts;

    for i = 1:nTotal
        numFrames = paramFrameCounts(i);
        frames = (startFrame:startFrame + numFrames - 1)';

        switch paramType{i}
            case 'Circular'
                theta = linspace(0, 2*pi, numFrames)';
                x = centerX + paramRadii(i) * cos(theta);
                y = centerY + paramRadii(i) * sin(theta);
                trajectoryKappaTrue{i} = (1/paramRadii(i)) * ones(numFrames,1);

            case 'GentleArc'
                arcSpanRad = deg2rad(paramArcSpanDeg(i));
                theta = linspace(0, arcSpanRad, numFrames)';
                x = centerX + paramRadii(i) * cos(theta);
                y = centerY + paramRadii(i) * sin(theta);
                trajectoryKappaTrue{i} = (1/paramRadii(i)) * ones(numFrames,1);

            case 'Elliptical'
                theta = linspace(0, 2*pi, numFrames)';
                x = centerX + paramSemiAxes(i,1) * cos(theta);
                y = centerY + paramSemiAxes(i,2) * sin(theta);
                a = paramSemiAxes(i,1); b = paramSemiAxes(i,2);
                trajectoryKappaTrue{i} = (a*b) ./ (a^2*sin(theta).^2 + b^2*cos(theta).^2).^1.5;

            case 'SpatiallyVarying'
                tt = linspace(0, paramSV.tSpan(i), numFrames)';
                f1 = paramSV.freq1(i); f2 = paramSV.freq2(i);
                a1 = paramSV.amp1(i);  a2 = paramSV.amp2(i);

                x = centerX + tt + a1*sin(f1*tt) + a2*sin(f2*tt);
                y = centerY + a1*1.2*sin(0.5*f1*tt) + a2*sin(1.1*f2*tt);

                dx  =  1 + a1*f1*cos(f1*tt) + a2*f2*cos(f2*tt);
                dy  =  a1*1.2*0.5*f1*cos(0.5*f1*tt) + a2*1.1*f2*cos(1.1*f2*tt);
                ddx = -a1*f1^2*sin(f1*tt) - a2*f2^2*sin(f2*tt);
                ddy = -a1*1.2*(0.5*f1)^2*sin(0.5*f1*tt) - a2*(1.1*f2)^2*sin(1.1*f2*tt);

                speedSq = dx.^2 + dy.^2;
                if min(speedSq) < 0.05^2
                    warning('Track %d (SpatiallyVarying): near-cusp detected (min speed = %.4f). Curvature may spike.', ...
                        i, sqrt(min(speedSq)));
                end

                trajectoryKappaTrue{i} = (dx.*ddy - dy.*ddx) ./ speedSq.^1.5;
                theta = tt;

            case 'FastVarying'
                tt = linspace(0, 1, numFrames)';   % normalized parameter
                nCycles = paramFV.nCycles(i);
                kAmp    = paramFV.kAmp(i);
                kBase   = paramFV.kBase(i);

                kappa_t = kBase + kAmp*sin(2*pi*nCycles*tt);   % prescribed curvature

                % Analytic integral of kappa_t (heading(0) = 0)
                heading = kBase*tt - (kAmp/(2*pi*nCycles))*cos(2*pi*nCycles*tt) ...
                          + (kAmp/(2*pi*nCycles));

                dx = cos(heading);
                dy = sin(heading);

                x = centerX + cumtrapz(tt, dx);
                y = centerY + cumtrapz(tt, dy);

                trajectoryKappaTrue{i} = kappa_t;
                theta = tt;
        end

        trajectoryTheta{i} = theta;

        if jitterStd > 0
            x = x + jitterStd * randn(numFrames, 1);
            y = y + jitterStd * randn(numFrames, 1);
        end
        if quantStep > 0
            x = round(x / quantStep) * quantStep;
            y = round(y / quantStep) * quantStep;
        end

        z = zValue * ones(numFrames, 1);
        trajectories{i} = [frames, x, y, z];
    end
end