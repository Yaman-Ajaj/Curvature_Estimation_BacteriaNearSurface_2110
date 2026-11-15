clear; clc;
raw = load('Ref0722_unfiltered.mat');   % <-- your unfiltered file
tracks = raw.AllBugs_unfilt;                  % adjust if the variable name differs

xs = []; 
for k = 1:min(numel(tracks),2000)
    M = tracks{k};
    if isempty(M), continue; end
    xs = [xs; M(~isnan(M(:,2)),2)]; %#ok<AGROW>
end
dx = abs(diff(sort(unique(round(xs,6)))));
fprintf('Smallest spacing between distinct x values: %.5f\n', min(dx));
fprintf('Median spacing: %.5f\n', median(dx));
fprintf('Distinct values / total samples: %d / %d\n', numel(unique(round(xs,6))), numel(xs));

for d = 2:5
    c = nchoosek(2*d, d);
    pooled = [];
    for k = 1:numel(tracks)
        M = tracks{k};
        if isempty(M), continue; end
        M = M(~any(isnan(M(:,2:3)),2),:);
        if size(M,1) < d+5, continue; end
        seg = cumsum([1; diff(M(:,1)) > 1]);
        for s = unique(seg)'
            m = M(seg==s, 2:3);
            if size(m,1) < d+3, continue; end
            pooled = [pooled; reshape(diff(m,d,1),[],1)]; %#ok<AGROW>
        end
    end
    lo = prctile(pooled,1); hi = prctile(pooled,99);
    tr = pooled(pooled>=lo & pooled<=hi);
    fprintf('d=%d  trimmed sigma = %.4f   plain sigma = %.4f\n', d, std(tr)/sqrt(c)/0.9342, std(pooled)/sqrt(c));

    sigTrim = std(tr)/sqrt(c)/0.9342;          % same correction as in your print line
    sdFit   = sigTrim*sqrt(c);                 % implied std of the d-th differences
    mu      = median(pooled);

    figure('Color','w'); hold on;
    view1 = pooled(pooled > prctile(pooled,0.1) & pooled < prctile(pooled,99.9));  % display range only
    histogram(view1, 'BinWidth', 0.325, 'Normalization','pdf', 'FaceColor',[0.6 0.6 0.6], 'EdgeColor','none');
    histogram(tr,    'BinWidth', 0.325, 'Normalization','pdf', 'FaceColor',[0.2 0.5 0.8], 'FaceAlpha',0.6, 'EdgeColor','none');

    xg = linspace(min(view1), max(view1), 400);
    plot(xg, normpdf(xg, mu, sdFit), 'r-', 'LineWidth', 1.5);
    xline(lo, 'k--'); xline(hi, 'k--');

    xlabel(sprintf('%d-th difference of position', d)); ylabel('probability density');
    title(sprintf('d = %d: trimmed \\sigma = %.4f (plain = %.4f)', d, sigTrim, std(pooled)/sqrt(c)));
    legend('all differences (0.1 to 99.9 pct)', 'kept after 1%/99% trim', 'Gaussian from trimmed \sigma', 'trim limits');
    grid on;
end
