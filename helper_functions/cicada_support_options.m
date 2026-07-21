function options = cicada_support_options(user_options)
% CICADA_SUPPORT_OPTIONS Return validated defaults for support/reliability products.
%
% options = cicada_support_options()
% options = cicada_support_options(user_options)
%
% The defaults preserve all historical CICADA outputs and add new products.
% No new product is used to alter denoising or IC classification.

options = struct();
options.Enabled = true;
options.Overwrite = false;
options.FailOnError = false;

% Component-derived evidence products.
options.SignalCoreProbabilities = [0.50 0.95];
options.NoiseCoreProbabilities = [0.50 0.95];
options.ExtendedSmoothingFWHMmm = 6;
options.ExtendedThresholdPercentRobustRange = 50;
options.WriteFixedThresholdExtendedMask = true;
options.FixedThresholdExtendedProbability = 0.50;

% Historical output provenance only. Historical commands remain in the
% original Auto/Manual scripts and are not replaced by these settings.
options.LegacySmoothingSigmaMM = 6;
options.LegacyThresholdPercentRobustRange = 50;

% Direct-data reliability products.
options.WriteStandardTSNR = true;
options.UseLowRobustTSNROutlierExclusion = true;
options.MinimumRobustTSNR = []; % no universal hard floor by default
options.RobustSigmaScale = 1.4826; % Gaussian consistency factor for MAD

if nargin < 1 || isempty(user_options)
    return
end
if ~isstruct(user_options)
    error('CICADA:SupportOptionsType', 'support_options must be a struct.');
end

fields = fieldnames(user_options);
known_fields = fieldnames(options);
unknown_fields = setdiff(fields, known_fields);
if ~isempty(unknown_fields)
    error('CICADA:UnknownSupportOption', 'Unknown support option(s): %s', ...
        strjoin(string(unknown_fields), ', '));
end
for i = 1:numel(fields)
    options.(fields{i}) = user_options.(fields{i});
end

validateattributes(options.Enabled, {'logical','numeric'}, {'scalar'});
validateattributes(options.Overwrite, {'logical','numeric'}, {'scalar'});
validateattributes(options.FailOnError, {'logical','numeric'}, {'scalar'});
validateattributes(options.SignalCoreProbabilities, {'numeric'}, {'vector','>=',0,'<=',1});
validateattributes(options.NoiseCoreProbabilities, {'numeric'}, {'vector','>=',0,'<=',1});
validateattributes(options.ExtendedSmoothingFWHMmm, {'numeric'}, {'scalar','nonnegative','finite'});
validateattributes(options.ExtendedThresholdPercentRobustRange, {'numeric'}, {'scalar','>',0,'<',100});
validateattributes(options.FixedThresholdExtendedProbability, {'numeric'}, {'scalar','>=',0,'<=',1});
validateattributes(options.WriteFixedThresholdExtendedMask, {'logical','numeric'}, {'scalar'});
validateattributes(options.WriteStandardTSNR, {'logical','numeric'}, {'scalar'});
validateattributes(options.UseLowRobustTSNROutlierExclusion, {'logical','numeric'}, {'scalar'});
validateattributes(options.RobustSigmaScale, {'numeric'}, {'scalar','positive','finite'});
if ~isempty(options.MinimumRobustTSNR)
    validateattributes(options.MinimumRobustTSNR, {'numeric'}, {'scalar','nonnegative','finite'});
end

options.Enabled = logical(options.Enabled);
options.Overwrite = logical(options.Overwrite);
options.FailOnError = logical(options.FailOnError);
options.WriteFixedThresholdExtendedMask = logical(options.WriteFixedThresholdExtendedMask);
options.WriteStandardTSNR = logical(options.WriteStandardTSNR);
options.UseLowRobustTSNROutlierExclusion = logical(options.UseLowRobustTSNROutlierExclusion);
options.SignalCoreProbabilities = unique(sort(options.SignalCoreProbabilities(:)'));
options.NoiseCoreProbabilities = unique(sort(options.NoiseCoreProbabilities(:)'));
if ~all(ismember([0.50 0.95], options.SignalCoreProbabilities)) || ...
        ~all(ismember([0.50 0.95], options.NoiseCoreProbabilities))
    error('CICADA:RequiredCoreProbabilities', ...
        'The rc1 group-QC pathway requires both 0.50 and 0.95 signal/noise core thresholds.');
end
end
