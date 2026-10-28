# Input and output conventions

## Catalog

Supply a FITS table in HDU 1. The existing selection is
`b_mode_mask == True` and `hsc_y3_zbin == 4`. The name `b_mode_mask` is a
catalog-selection flag; this workflow does not compute a B-mode null test.

| Column | Meaning |
| --- | --- |
| `i_ra`, `i_dec` | RA and declination in degrees |
| `i_hsmshaperegauss_derived_weight` | Galaxy weight |
| `b_mode_mask`, `hsc_y3_zbin` | Selection fields |
| Ellipticity components | Required for noise calibration; see below |

For each component `a=1,2`, the noise loader chooses the first existing column
from `i_hsmshaperegauss_e{a}`, `i_hsmshaperegauss_derived_e{a}`,
`i_hsmshaperegauss_shape_{a}`, `e{a}`. It rejects nonfinite ellipticities for
otherwise selected valid galaxies. The simulation/noise loader filters
nonfinite positions/weights and nonpositive weights; the preserved metadata
mapper uses the selected catalog directly. Input quality must therefore be
consistent across these stages. The tested XMM sample contains 300101 galaxies.

Noise is calibrated from weighted, uncentered ellipticity second moments:
`R_a = 1 - <e_a^2>_w` and
`sigma_component_squared = (<(e_1/(2R_1))^2>_w + <(e_2/(2R_2))^2>_w)/2`.
This is the variance of **one** shear component, not the sum over components.
The calibration routine is shared with the noise simulation script, but
importing it for analytic covariance does not run simulations.

## Spectrum and angular bins

- Spectrum NPZ: a one-dimensional `cl_theory` array, or a key explicitly
  supplied using `--cl-key`. Array index is ell, including ell=0 and 1.
  Supply at least 1536 entries. The model uses entries through ell=1535,
  sets ell=0,1 to zero, and treats the rest directly as EE. BB is zero.
- The input is `C_ell`, not `ell(ell+1)C_ell/(2 pi)`. If its source is a
  convergence calculation, confirm that directly reusing its values as EE is
  the intended convention; the program performs no conversion.
- Bin-settings NPZ: `theta` (length 10) and `theta_edges` (length 11), in
  **degrees**. Supply positive, increasing, logarithmically spaced edges for
  the preserved TreeCorr setup. The metadata code reads `theta` but constructs
  centers as `sqrt(theta_edges[:-1] * theta_edges[1:])`.
- Internal spin-2 transforms use radians. Saved angular-bin values use degrees;
  comparison plots display arcminutes by multiplying by 60.

The spectrum and bin-settings archives are user-supplied inputs; they need not
contain simulation data. They can be written with `np.savez` from the user's
chosen arrays. No cosmological-parameter or redshift-distribution calculation
is performed by this package.

## Saved products

| Product / key | Meaning |
| --- | --- |
| `pair_meta.npz`, merged pair NPY files | Pixel weights/geometry and unordered pixel pairs per angular bin |
| `model.npz: prediction_mean` | Twenty predicted, weighted grid-estimator correlations (all plus, then minus) |
| `model.npz: xi_plus, xi_minus` | Correlation functions on `theta_grid` in radians; these are not ten-bin means |
| `signal_ht.npz: covariance_joint` | Signal-only covariance, shape (20,20) |
| `noise_ht.npz: covariance_joint` | Total signal + shape-noise covariance, shape (20,20) |
| `covariance_sem` | Standard error of the mean across the HT thinning realizations |
| `covariance_realizations` | Individual HT estimates, shape (100,20,20) for the default run |
| `covariance_noise_added` | Total minus matched signal-only HT covariance: mixed term **plus** noise-only term |
| `pixel_noise_variance` | Weighted-mean noise variance per component per active pixel |
| `full_weight_sums`, `retained_counts`, `bin_seeds` | Normalization and thinning provenance |

`noise_ht.npz` is a suggested output name; it contains the total covariance,
not a noise-only covariance. Use `--signal-only-reference` with the matching
signal-only HT file to fill `covariance_noise_added`. Without a matching
baseline it may be empty. Without `--sims`, the stored `simulation_covariance_*`
arrays are empty and `simulation_nsims` is zero; this means no simulation
reference was supplied, not that the physical covariance is zero.

For each matrix, indices 0..9 represent plus bins and 10..19 minus bins.
The cross blocks satisfy `C_minus_plus = C_plus_minus.T`; a cross block need
not be symmetric on its own. `offset=0` compares equal angular-bin indices
within each block, while `offset=1` compares neighboring angular bins.

## File locations and restarts

Use one fresh `WORK` directory per configuration. Prepared models store
absolute pair-file paths and hashes; keep the model and pairs at those paths
while running the documented workflow. Regenerate the model if relocating the
run for the noise workflow. Retain the directory structure of the code because
neighboring folders share modules. Final covariance/simulation files are not
overwritten; matching checkpoints support resuming. Start a new output path
when changing scientific inputs/settings. The preprocessing scripts should
also be run into a fresh directory, as instructed in the quickstart.
