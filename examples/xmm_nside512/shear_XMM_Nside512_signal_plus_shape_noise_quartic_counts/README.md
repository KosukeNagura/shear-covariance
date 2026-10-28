# Shear Nside512 grid covariance with shape noise

This implements the closed HEALPix-grid estimator for signal plus independent catalog shape noise. It reuses the signal-only model, pair lists, pixel weights, and covariance kernel. For each active pixel i, the catalog-calibrated per-component noise variance is

`v_i = sigma_component_squared * sum_g(w_g^2) / sum_g(w_g)^2`.

The code adds `v_i I_2` to the coincident-pixel shear covariance block `S_ii`. The existing four-block Wick kernel then includes signal-only, signal-noise, and noise-only contractions together. This is the grid estimator; TreeCorr catalog simulations are optional comparison products, and the estimators must not be conflated.

Outputs refuse overwrite and support checkpoints. Pass `--signal-only-reference /path/to/signal_ht.npz` to populate `covariance_noise_added`; the model hash and thinning seeds/counts must match exactly. Without this option, the legacy sibling-output lookup is retained. For HT this lets the paired difference isolate noise contributions with lower Monte Carlo scatter.

Run-time inputs (kept outside the source package):
- signal-only `outputs/model.npz` and verified pair files;
- the authorized XMM catalog including ellipticity columns.

By default, the code calibrates the noise variance directly from the catalog by reusing the established simulation loader and `noise_parameters` function. No simulation is executed or required. For compatibility, `--sims /path/to/simulation.npz` retains the prior archive-input route. The two routes use the same variance convention. In direct-catalog mode, `simulation_covariance_*` arrays are empty, `simulation_reference` is empty, and `simulation_nsims=0`; they must not be interpreted as simulated zero covariances.

For the current validation, run only the 100-realization HT estimate with keep probability 0.001 from the down-sampling folder. The non-thinned/full-row calculation and row merger are outside this workflow. Site-specific validation PBS examples are under `validation/`; adapt resource requests for another system. Check the inputs and output location before submitting, and do not mix files from different configurations.
