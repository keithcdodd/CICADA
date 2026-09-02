function tag = cicada_smoothing_fwhm_tag(fwhm_mm)
%CICADA_SMOOTHING_FWHM_TAG Make a filename-safe, non-rounded FWHM tag.
%   Integer values retain the historical form (for example, s6). Decimal
%   points and exponent signs are encoded so distinct non-integer requests
%   are not silently assigned the same rounded output name.

validateattributes(fwhm_mm, {'numeric'}, ...
    {'real','scalar','finite','nonnegative'}, mfilename, 'fwhm_mm');

value_text = sprintf('%.15g', double(fwhm_mm));
value_text = strrep(value_text, '.', 'p');
value_text = strrep(value_text, '+', '');
value_text = strrep(value_text, '-', 'm');
tag = ['s', value_text];
end
