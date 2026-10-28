# XMM Nside=512 example: current implementation

This is the current shear two-point covariance workflow for the XMM field, one of the regions in the HSC Y3 catalog used in this analysis. It uses the bin-4 sample and Nside=512. It reuses pair products made by the earlier catalog precomputation and includes the original Nside=512 metadata, pair-worker, and merge scripts.

## Pipeline

1. Generate `pair_meta.npz` and `pairs_merged/pairs_bin_XX.npy` from the catalog using the included established metadata/pair scripts.
2. `shear_XMM_Nside512_signal_only_quartic_counts/prepare_shear.py` builds the shear model from the input spectrum and those pair lists.
3. Optional comparison simulations: `shear_XMM_Nside512_signal_only_sims/simulate_shear.py` for signal-only grid/catalog measurements, and `shear_XMM_Nside512_signal_plus_shape_noise_sims/simulate_shear_shape_noise.py` for signal/noise/total catalog measurements. Neither is required to compute analytic covariance.
4. `shear_XMM_Nside512_signal_plus_shape_noise_quartic_counts/compute_covariance_noise.py` contains the existing covariance implementation; the validation workflow here runs only the 100-realization thinned estimate (keep probability 0.001). The non-thinned calculation is not part of this workflow.
5. The notebook under `shear_XMM_Nside512_signal_plus_shape_noise_quartic_counts_down_sampling/notebooks/` compares covariance blocks and the catalog simulation.

## Existing pair-preparation products

The current workflow already used catalog-derived Nside=512 pair products from the earlier XMM precomputation: `pair_meta.npz` and `pairs_merged/pairs_bin_XX.npy`. Regenerate those products with the included original implementation for a complete reproduction; historical products serve as comparison references. Keep them in the `outputs/precompute/` directory expected by `prepare_shear.py`.

The original Nside=512 metadata, pair-worker, and merge scripts from `local/XMM/Nside512/` are included under `pair_precompute_existing/`; the original pair-chunk submission wrapper is adapted there as well. The pair-worker and merge algorithms are preserved; the config adapter replaces machine-specific paths with environment variables. The required original `galaxy_catalog_to_healpix` utility is included in `legacy_2pcf_funcs/`. No external `2PCF/funcs` checkout is required.

To build the model from the existing precomputation, run from this example directory:

```sh
python shear_XMM_Nside512_signal_only_quartic_counts/prepare_shear.py \
  --pair-base shear_XMM_Nside512_signal_only_quartic_counts/outputs/precompute \
  --cl-file /path/to/authorized/spectrum.npz --cl-key cl_theory
```

## Remaining study-specific inputs and limitations

The signal spectrum is an external calculation input. Historical model/simulation results are only comparison references. Simulation checks against prior signal-only/convergence results are optional and explicitly requested with reference flags. No catalog, spectrum, or generated NPZ pair products belong in Git unless their data-use terms explicitly permit redistribution.

The catalog simulation measures ξ± directly with TreeCorr. The analytic covariance code evaluates the weighted HEALPix-grid pair estimator. These are distinct estimators and the notebook presents their comparison as a diagnostic. This snapshot preserves the current study setup; it is not yet a fully survey-independent package.
