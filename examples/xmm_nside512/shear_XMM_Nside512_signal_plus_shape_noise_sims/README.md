# Catalog simulation with shape noise

This reproduces the current shear two-point catalog simulation: EE signal maps, independent isotropic Gaussian noise on each galaxy shear component, and TreeCorr measurements of signal, noise, and total shear. It is a simulation validation companion to the grid analytic covariance code.

## Scientific choices preserved from the current XMM run

- XMM selection in tomographic bin 4, Nside=512, and the ten angular bins in the model file.
- EE spectrum from `model.npz`, BB=0, ell through 1535, `pixwin=False`, and no smoothing.
- The catalog-based measurement uses the TreeCorr default metric/slop settings.
- Shape noise is independent per galaxy and per component. Its variance is calibrated using weighted ellipticity moments and the responsivities used in the convergence analysis; an optional regression check compares it with an explicitly supplied convergence reference.
- When `--reference-signal` is supplied, each signal realization is checked against that catalog realization for the same seed.

The code intentionally does not change the estimators or signal convention from the current analysis. For a different survey, first verify the catalog schema, cuts, weighting, noise calibration, and model correspondence.

## Inputs

Pass authorized files with CLI options; no catalog is bundled with the package. The calculation needs the model (`--model`) and catalog (`--catalog`). Prior signal-only and convergence results are optional validation inputs (`--reference-signal`, `--reference-convergence`); they are never used to generate signal or noise. Omitting them records `NOT_REQUESTED` in the log and null reference hashes, rather than claiming a passed reference check. The required columns and catalog selection are documented in `simulate_shear_shape_noise.py`.

## Running

Activate an environment with the versions in `requirements-current.txt`, then use the small pilot first:

```sh
python simulate_shear_shape_noise.py --model /path/to/new_analysis_run/model.npz --catalog /path/to/authorized/XMM.fits --nsims 3 --threads 8 --output /path/to/new_analysis_run/noise_3sims.npz
```

Review the saved configuration and, if references were supplied, the reference-check messages. The command above uses no historical reference. Add either optional reference flag only when performing the corresponding historical check. The numerical simulation is identical. The integration test now compares every saved scientific array between runs with and without references. The idark reproduction example in `validation/run_from_catalog.pbs` uses no historical simulation inputs and compares outputs afterward.

The output contains per-realization `xi_signal_sims`, `xi_noise_sims`, and `xi_total_sims`, their means and sample covariance matrices, imaginary components, seed and input hashes, noise calibration, and configuration.
