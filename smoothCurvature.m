function kSmoothed = smoothCurvature(k, windowSize, method)
if nargin < 2
 windowSize = 5;
end
if nargin < 3
 method = 'movmean';
end

if windowSize <= 1
    kSmoothed = k;
    return;
end

n = numel(k);

% Guard: filters requiring a fixed frame length can't operate when the
% input is shorter than the requested window. Shrink the window to fit
% (kept odd for sgolay), or fall back to no smoothing if too short even
% for a minimal window.
if strcmpi(method, 'sgolay') || strcmpi(method, 'median')
    if n < 5
        kSmoothed = k; % too short to smooth meaningfully at all
        return;
    end
    if windowSize > n
        windowSize = n;
        if mod(windowSize, 2) == 0
            windowSize = windowSize - 1; % keep odd, still fits within n
        end
    end
end


switch lower(method)
case 'movmean'
 kSmoothed = movmean(k, windowSize, 'omitnan');

case 'median'
% medfilt1 doesn't have a built-in 'omitnan' option, so NaNs are
% temporarily filled via linear interpolation for the filtering step
% only, then re-masked back to NaN at their original positions
% afterward, so gaps aren't silently invented.
 nanMask = isnan(k);
if any(nanMask) && ~all(nanMask)
 kFilled = k;
 kFilled(nanMask) = interp1(find(~nanMask), k(~nanMask), find(nanMask), 'linear', 'extrap');
else
 kFilled = k;
end
 kSmoothed = medfilt1(kFilled, windowSize);
 kSmoothed(nanMask) = NaN;

case 'sgolay'
 polyOrder = 3;
if mod(windowSize, 2) == 0
 windowSize = windowSize + 1; % sgolayfilt requires an odd frame length
end
if windowSize <= polyOrder
 windowSize = polyOrder + 2; % must exceed polynomial order; keep it odd
if mod(windowSize, 2) == 0
 windowSize = windowSize + 1;
end
end

% Same NaN-fill/re-mask approach as the median case, since sgolayfilt
% also has no native NaN handling
 nanMask = isnan(k);
if any(nanMask) && ~all(nanMask)
 kFilled = k;
 kFilled(nanMask) = interp1(find(~nanMask), k(~nanMask), find(nanMask), 'linear', 'extrap');
else
 kFilled = k;
end
 kSmoothed = sgolayfilt(kFilled, polyOrder, windowSize);
 kSmoothed(nanMask) = NaN;

otherwise
 error('smoothCurvature:unknownMethod', ...
 'Unknown smoothing method "%s". Use ''movmean'', ''median'', or ''sgolay''.', method);
end
end