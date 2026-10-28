# Signal-only validation simulation

`simulate_shear.py` generates Gaussian EE shear maps, measures the configured grid estimator, and (unless `--grid-only` is set) measures the selected galaxy catalog with TreeCorr. It is a regression/estimator-comparison tool; the noise-inclusive counterpart is in `shear_XMM_Nside512_signal_plus_shape_noise_sims/`.

The run needs `outputs/model.npz` from the common signal-only folder and catalog-derived pair lists. Pass the authorized catalog path with `--catalog`. This is the original study implementation and retains its model conventions, nearest-pixel shear assignment, and TreeCorr defaults. It imports the existing `load_pairs` helper from the common signal-only covariance module; keep the sibling-folder layout intact.

Cluster submission files from the original idark environment are intentionally not included. Adapt a scheduler script to the target cluster, resource limits, environment activation, and local paths. Do not submit a large run before a small pilot passes.
