% synthCurvatureSpeedTest.m
% Tests the curvature-vs-speed power-law fit on synthetic circling tracks with known
% truth. Scenario A: curvature fixed, speed varied (true exponent b = 0), which measures
% the artifact that noise alone creates. Scenario B: circling rate omega fixed, so
% k = omega/v (true b = -1), which checks whether the pipeline recovers it. Tracks get
% the same noise as the real data (Gaussian jitter + pixel quantisation), the same
% spline-smoothed tangent-angle estimator and the same parameters as analyzeAllTracks.
% the script generates circling tracks with a known exponent, using real speeds and 
% section lengths, and compares the fitted exponent with the true one.

clear; 
%clc;

%% Settings
nTracks     = 1500;            % synthetic tracks per repetition
nRep        = 5;               % repetitions with different seeds
omega0      = 0.83/15;         % median circling rate, rad per frame (real median / fps)
sigma = 0.25;  tolMult = 1.29;  step = 7;     % same as analyzeAllTracks
quantStep   = 0.325;
jitter      = sqrt(sigma^2 - quantStep^2/12);
names       = {'A','B','C','D'};
bTrueList   = [0 -1 -0.7 -1.3];                % true exponent of each scenario
kSpreadList = [0.5 1.05];                       % lognormal spread of k at fixed speed

Rr   = load('realSections.mat');               % real speeds and lengths (curvatureVsSpeed.m)
vMed = median(Rr.vReal);

%% Run
for sp = kSpreadList
    for sc = 1:numel(bTrueList)
        bTrue = bTrueList(sc);
        bEst = nan(nRep,1);  bIdeal = nan(nRep,1);  ratio = nan(nRep,1);

        for rep = 1:nRep
            rng(100*rep + sc);
            j  = randi(numel(Rr.vReal), nTracks, 1);       % real (speed, length) pairs
            v  = Rr.vReal(j);   nF = round(Rr.nFrReal(j));

            % true curvature: power law in speed times scatter; omega = k*v
            kTrue = (omega0/vMed) * (v/vMed).^bTrue .* exp(sp*randn(nTracks,1));
            w     = kTrue .* v;                            % rad per frame

            kMean = nan(nTracks,1);
            for i = 1:nTracks
                if nF(i) < 2*step + 3, continue; end
                th = w(i) * (0:nF(i)-1)';
                x  = cumsum(v(i)*cos(th));   y = cumsum(v(i)*sin(th));
                x  = round((x + jitter*randn(nF(i),1)) / quantStep) * quantStep;
                y  = round((y + jitter*randn(nF(i),1)) / quantStep) * quantStep;
                [k, ~] = curvature_estimation_splineSmoothed([x y], sigma, step, (1:nF(i))', tolMult);
                k = k(isfinite(k));
                if ~isempty(k), kMean(i) = mean(abs(k)); end
            end

            ok = isfinite(kMean) & kMean > 0;
            pE = polyfit(log(v(ok)), log(kMean(ok)), 1);   % fit on estimated curvature
            pI = polyfit(log(v(ok)), log(kTrue(ok)), 1);   % same fit on TRUE curvature
            bEst(rep) = pE(1);  bIdeal(rep) = pI(1);
            ratio(rep) = median(kMean(ok) ./ kTrue(ok));
            if rep == 1
                r  = kMean(ok) ./ kTrue(ok);   vv = v(ok);
                fprintf('   est/true k percentiles 5/25/50/75/95: %.2f %.2f %.2f %.2f %.2f\n', prctile(r,[5 25 50 75 95]));
                q  = prctile(vv, [0 33 67 100]);   [~,~,bn] = histcounts(vv, q);
                fprintf('   median est/true k by speed (slow, mid, fast): %.2f %.2f %.2f\n', ...
                    median(r(bn==1)), median(r(bn==2)), median(r(bn==3)));
                pR = robustfit(log(vv), log(kMean(ok)));
                fprintf('   robust fit b = %+.2f\n', pR(2));
            end
        end

        fprintf('spread %.1f | %s true b %+.1f | ideal fit %+.2f | estimated fit %+.2f (sd %.2f) | bias %+.2f | est/true k %.2f\n', ...
            sp, names{sc}, bTrue, mean(bIdeal), mean(bEst), std(bEst), mean(bEst - bIdeal), mean(ratio));
    end
    fprintf('\n');
end