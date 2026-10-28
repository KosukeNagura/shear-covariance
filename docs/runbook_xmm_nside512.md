# Optional simulations and comparison plots

First complete the [analytic quickstart](analytic_quickstart.md). It defines
`EXAMPLE`, `CATALOG` and `WORK`, and creates the model and both HT outputs.
Run the following on a compute allocation with the same environment. No
historical results are needed for generation.

## Pilot and production simulations

Run a short pilot with separate filenames:

```bash
python "$EXAMPLE/shear_XMM_Nside512_signal_only_sims/simulate_shear.py" \
  --model "$WORK/model.npz" --catalog "$CATALOG" \
  --nsims 3 --threads 8 --output "$WORK/signal_3sims.npz"

python "$EXAMPLE/shear_XMM_Nside512_signal_plus_shape_noise_sims/simulate_shear_shape_noise.py" \
  --model "$WORK/model.npz" --catalog "$CATALOG" \
  --nsims 3 --threads 8 --output "$WORK/noise_3sims.npz"
```

After the pilot completes successfully, run production with new filenames:

```bash
python "$EXAMPLE/shear_XMM_Nside512_signal_only_sims/simulate_shear.py" \
  --model "$WORK/model.npz" --catalog "$CATALOG" \
  --nsims 3000 --threads 8 --output "$WORK/signal_3000sims.npz"

python "$EXAMPLE/shear_XMM_Nside512_signal_plus_shape_noise_sims/simulate_shear_shape_noise.py" \
  --model "$WORK/model.npz" --catalog "$CATALOG" \
  --nsims 3000 --threads 8 --output "$WORK/noise_3000sims.npz"
```

Request at least the CPU count passed to `--threads`. Both commands generate
their own signal maps using the established seeds. The noise command records
signal-only, noise-only and total catalog measurements; its total is measured
from signal plus noise at each galaxy, including mixed contributions. The
signal-only command also records grid measurements. The optional
`--reference-signal` and `--reference-convergence` flags belong only to the
noise script and enable explicit historical checks, not input generation.
`--require-conda` is an idark convenience guard for the environment name
`cellcov_3`; it is not required for other compatible environments.

## Plot the shape-noise comparison

Open the notebook under
`examples/xmm_nside512/shear_XMM_Nside512_signal_plus_shape_noise_quartic_counts_down_sampling/notebooks/compare_shear_shape_noise_ht_with_catalog_sims_nside512.ipynb`.
In its first code cell, set both absolute paths:

```python
SIM_FILE = "/path/to/new_analysis_run/noise_3000sims.npz"
HT_FILE = "/path/to/new_analysis_run/noise_ht.npz"
```

Restart the kernel and run all cells. The notebook requires 3000 simulations,
100 HT realizations, p=0.001 and matching model hashes. It prints the input
paths. The historical automatic path lookup remains available when both path
variables are None; explicit paths avoid selecting old results accidentally.

Figures show four covariance blocks at offsets 0 and 1, correlation matrices,
and the matched total-minus-signal covariance increment. The blue error bars
use a Gaussian approximation to sample-covariance uncertainty; orange error
bars show HT mean uncertainty. They quantify different sources of error.
Catalog and grid estimators differ as described in [the method](method.md).
The validation notebook under `validation/notebooks/` has the original idark
job paths prefilled and is intended only for that reproduction record.
