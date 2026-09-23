## CICADA v1.2.0

## Group CICADA processing updates
- Group CICADA now uses matched temporal processing by default: the same detrending and optional temporal-filtering transform is applied to the original BOLD data and complete MELODIC mixing matrix before nonaggressive ICA cleanup.
- Historical sequential processing remains available with temporal_cleanup_mode = 'legacy_sequential'.
- Mask-normalized Gaussian smoothing is now used at functional-mask boundaries to avoid attenuation from outside-mask zeros.
- Temporal filtering validation, cutoff/Nyquist/run-duration checks, and processing provenance have been strengthened.
- Degree-2 polynomial detrending remains the default.
- Temporal bandpass filtering remains optional and is off by default.
- Automatic and Manual CICADA IC-classification behavior is unchanged.

## Documentation
- Updated Group CICADA methods, parameter documentation, and manuscript-language examples to reflect the current processing architecture.

These changes were regression-tested against prior CICADA outputs, including exact reproduction of previously validated Group CICADA prepared images where behavior was intended to remain unchanged.