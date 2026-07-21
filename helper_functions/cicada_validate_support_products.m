function summary = cicada_validate_support_products(task_dir, cicada_type)
% CICADA_VALIDATE_SUPPORT_PRODUCTS Validate nesting, ranges, and geometry.
products = cicada_locate_support_products(task_dir, cicada_type);
required = {'FuncMaskConstrained','DirectReliabilityMask','RobustTSNR', ...
    'SignalICOverlapRaw','NoiseICOverlapRaw','SignalCoreP50','SignalCoreP95', ...
    'SignalExtendedAdaptive','SignalNoiseBalance','AnyICEvidence'};
missing = strings(0,1);
for i = 1:numel(required)
    if ~isfield(products, required{i}) || isempty(products.(required{i})) || ~isfile(products.(required{i}))
        missing(end+1,1) = string(required{i}); %#ok<AGROW>
    end
end
if ~isempty(missing)
    error('CICADA:MissingSupportProducts', 'Missing support products: %s', strjoin(missing, ', '));
end

geometry_paths = cellfun(@(f) products.(f), required, 'UniformOutput', false);
cicada_assert_nifti_geometry(products.FuncMaskConstrained, geometry_paths{:});

base = niftiread(products.FuncMaskConstrained) > 0;
direct = niftiread(products.DirectReliabilityMask) > 0;
rtsnr = single(niftiread(products.RobustTSNR));
signal = single(niftiread(products.SignalICOverlapRaw));
noise = single(niftiread(products.NoiseICOverlapRaw));
s50 = niftiread(products.SignalCoreP50) > 0;
s95 = niftiread(products.SignalCoreP95) > 0;
extended = niftiread(products.SignalExtendedAdaptive) > 0;
balance = single(niftiread(products.SignalNoiseBalance));
any_evidence = single(niftiread(products.AnyICEvidence));

arrays = {direct,rtsnr,signal,noise,s50,s95,extended,balance,any_evidence};
for i = 1:numel(arrays)
    if ~isequal(size(base), size(arrays{i}))
        error('CICADA:SupportGeometryMismatch', 'Support product %d has mismatched dimensions.', i);
    end
end
if any(direct(:) & ~base(:)) || any(s50(:) & ~base(:)) || ...
        any(s95(:) & ~base(:)) || any(extended(:) & ~base(:))
    error('CICADA:SupportMaskNotNested', 'A support mask extends outside funcmask_constrained.');
end
if any(s95(:) & ~s50(:))
    error('CICADA:CoreMaskNesting', 'The p95 signal core is not a subset of p50.');
end
if any(~isfinite(rtsnr(:))) || any(rtsnr(:) < 0)
    error('CICADA:InvalidRobustTSNR', 'Robust tSNR contains nonfinite or negative values.');
end
if any(~isfinite(signal(:))) || any(signal(:) < -1e-6 | signal(:) > 1+1e-6) || ...
        any(~isfinite(noise(:))) || any(noise(:) < -1e-6 | noise(:) > 1+1e-6)
    error('CICADA:InvalidComponentEvidence', 'Raw component evidence is outside [0,1] or nonfinite.');
end
if any(~isfinite(balance(:))) || any(balance(:) < -1-1e-6 | balance(:) > 1+1e-6) || ...
        any(~isfinite(any_evidence(:))) || any(any_evidence(:) < -1e-6 | any_evidence(:) > 1+1e-6)
    error('CICADA:InvalidDerivedComponentEvidence', ...
        'Derived component-evidence maps have invalid ranges or nonfinite values.');
end

summary = table(string(task_dir), string(cicada_type), nnz(base), nnz(direct), ...
    nnz(direct)/max(nnz(base),1), nnz(s50), nnz(s95), nnz(extended), ...
    median(double(rtsnr(base & rtsnr>0)), 'omitnan'), ...
    median(double(signal(base)), 'omitnan'), median(double(noise(base)), 'omitnan'), true, ...
    'VariableNames', {'TaskDirectory','CICADAType','ConstrainedVoxels', ...
    'DirectReliableVoxels','DirectReliableFraction','SignalCoreP50Voxels', ...
    'SignalCoreP95Voxels','SignalExtendedVoxels','MedianRobustTSNR', ...
    'MedianSignalICEvidence','MedianNoiseICEvidence','ValidationPassed'});
end
