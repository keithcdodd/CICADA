function [keep, rtsnr, valid, lowCut, fixedCut, cutoff] = ...
    cicada_robust_tsnr_selection(mask, tmedian, tmad, options)
% Select reliable voxels using CICADA's existing robust-tSNR rule.
% This helper performs no file I/O.

mask = logical(mask);

assert(isequal(size(mask), size(tmedian), size(tmad)), ...
    'Robust-tSNR inputs must have matching dimensions.');

robust_sigma = options.RobustSigmaScale .* tmad;

valid = mask & isfinite(tmedian) & ...
    isfinite(robust_sigma) & robust_sigma > 0;

rtsnr = zeros(size(tmedian), 'single');

rtsnr(valid) = ...
    tmedian(valid) ./ robust_sigma(valid);

rtsnr(~isfinite(rtsnr) | rtsnr < 0) = 0;

vals = double(rtsnr(valid & rtsnr > 0));

if isempty(vals)
    error('CICADA:NoValidRobustTSNR', ...
        'No valid robust-tSNR voxels in mask.');
end

lowCut = 0;

if options.UseLowRobustTSNROutlierExclusion
    log_vals = log1p(vals);

    low_flags = isoutlier(log_vals, 'median') & ...
        log_vals < median(log_vals, 'omitnan');

    if any(low_flags)
        lowCut = max(vals(low_flags));
    end
end

fixedCut = 0;

if ~isempty(options.MinimumRobustTSNR)
    fixedCut = options.MinimumRobustTSNR;
end

cutoff = max(lowCut, fixedCut);

keep = valid & rtsnr > cutoff;
end