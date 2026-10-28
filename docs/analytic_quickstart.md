# Analytic shear covariance from a catalog and a spectrum

This example preserves the validated XMM HSC Y3 bin-4, Nside512 conventions.
It is not yet a generic catalog-schema or arbitrary-tomography interface.
Use the original catalog column names and selection described in the code.
No simulation result is required by the analytic calculation.

Inputs are an authorized FITS catalog, an NPZ containing `cl_theory` (C_ell
indexed from ell=0), and an NPZ with the ten-bin `theta` / `theta_edges` settings
in degrees. `theta_edges` has length 11. The model retains ell through 1535,
EE only, with the original signal conventions. Use a fresh output directory.

Activate an environment with `requirements-current.txt` before running. The
commands below are for a compute allocation; do not run the heavy calculations
on a shared login node. Replace the example paths with actual absolute paths.

```bash
EXAMPLE='/path/to/shear-covariance-reproducibility/examples/xmm_nside512'
CATALOG='/path/to/XMM.fits'
SPECTRUM='/path/to/spectrum.npz'
BINS='/path/to/angular_bin_settings.npz'
WORK='/path/to/new_analysis_run'
mkdir "$WORK"

export SHEAR_PAIR_CATALOG="$CATALOG"
export SHEAR_PAIR_THETA_NPZ="$BINS"
export SHEAR_PAIR_WORK_ROOT="$WORK/pairs"
export SHEAR_2PCF_FUNCS="$EXAMPLE/pair_precompute_existing/legacy_2pcf_funcs"
cd "$EXAMPLE/pair_precompute_existing"
python prep_pair_metadata.py
```

Run all original pair workers before merging. This serial orchestration uses
the established worker without changing its algorithm:

```bash
python - <<'PY'
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
base = Path(os.environ['SHEAR_PAIR_WORK_ROOT'])/'outputs/precompute'
with np.load(base/'pair_meta.npz') as f:
    n = len(f['pix_idx'])
for chunk, start in enumerate(range(0, n, 2000)):
    subprocess.run([sys.executable, 'pair_chunk_worker.py', str(chunk),
                    str(start), str(min(start+2000, n))], check=True)
subprocess.run([sys.executable, 'merge_pair_chunks.py'], check=True)
PY

python "$EXAMPLE/shear_XMM_Nside512_signal_only_quartic_counts/prepare_shear.py" \
  --pair-base "$WORK/pairs/outputs/precompute" \
  --cl-file "$SPECTRUM" --cl-key cl_theory \
  --output "$WORK/model.npz"
```

Compute both thinned analytic covariances with the same sampling settings:

```bash
python "$EXAMPLE/shear_XMM_Nside512_signal_only_quartic_counts/compute_covariance.py" \
  --model "$WORK/model.npz" --mode thinned --realizations 100 \
  --base-seed 20260828 --left-chunk 0 --right-chunk 2000 \
  --output "$WORK/signal_ht.npz"

python "$EXAMPLE/shear_XMM_Nside512_signal_plus_shape_noise_quartic_counts/compute_covariance_noise.py" \
  --model "$WORK/model.npz" --catalog "$CATALOG" \
  --mode thinned --realizations 100 --base-seed 20260828 \
  --left-chunk 0 --right-chunk 2000 \
  --signal-only-reference "$WORK/signal_ht.npz" \
  --output "$WORK/noise_ht.npz"
```

The model contains predicted grid-estimator xi and its required intermediate
correlations. The two HT files contain 20x20 covariance matrices ordered as ten
plus bins then ten minus bins, their thinning SEM and per-realization results.
The noise-added file also stores the matched signal-only subtraction. It has no
simulation covariance unless an optional `--sims` archive was supplied.

The two simulation scripts and comparison notebook remain included for separate
validation. Do not label a plot made with existing simulation outputs as a fresh
simulation reproduction. Catalog and grid estimators differ as described in
`method.md`. No non-thinned production calculation is required by this workflow.
