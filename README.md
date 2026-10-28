# Shear covariance

Compute the Gaussian covariance of shear two-point correlation functions
(`xi_plus`, `xi_minus`), with or without shape noise, using weighted HEALPix
pixel pairs and pair thinning. Optional catalog simulations and comparison
notebooks are included.

This release preserves the implementation tested on the **XMM field of HSC Y3,
tomographic bin 4, Nside=512**. It is a collection of runnable scripts, not a
`pip install` package or a general interface for arbitrary survey catalogs.

## Inputs and outputs

| Input | Purpose |
| --- | --- |
| HSC-style FITS catalog | Galaxy positions, weights, selection and ellipticity-based noise calibration |
| Theoretical power spectrum (`C_ell`) in NPZ | EE signal input; array index is ell |
| Angular-bin settings in NPZ | Ten bin centers and eleven edges, in degrees |

The pipeline generates pixel metadata, pair lists, a signal model containing
predicted grid correlations, and a **20 x 20 covariance matrix** ordered as
all ten plus bins followed by all ten minus bins. Thinning uncertainty and
per-realization results are saved alongside the mean covariance.

Catalogs and generated products are not bundled. Historical simulation or
model outputs are **not required to run the analytic calculation**. They are
used only as regression references in the study-specific validation scripts.

## Start here

1. Read [input/output conventions](docs/inputs_and_outputs.md).
2. Follow [the analytic quickstart](docs/analytic_quickstart.md) to run
   catalog → pairs → model → thinned covariance, with and without noise.
3. Optionally run [catalog simulations and plots](docs/runbook_xmm_nside512.md).
4. Read [the method summary](docs/method.md) and
   [reproduction results](docs/reproducibility_validation.md).

Use an existing compatible scientific Python environment or prepare a separate
one using `requirements-current.txt`. It records selected tested package
versions; it is not a complete environment lockfile. No environment updates
are performed by these scripts. `pytest` is optional developer tooling and is
not needed for the analysis commands. Run the heavy stages on compute nodes.

## Scope and interpretation

- The documented workflow uses **thinned** covariance: pair retention
  probability 0.001, 100 realizations by default. Unthinned production runs are
  not required. Existing full-mode code is retained for implementation continuity.
- The input spectrum is treated directly as spin-2 EE, with BB=0,
  ell_max=1535, no pixel-window factor and no smoothing. No convergence-to-shear
  spectrum conversion is applied.
- Analytic results describe a **grid pair estimator**. The comparison
  simulations measure **catalog galaxy pairs with TreeCorr**. Their agreement
  is a diagnostic, not a claim that the two estimators are identical.
- The current selection, noise calibration and simulation guards target the
  XMM setup above. Other redshift bins, catalogs or B-mode statistics require
  explicit adaptation and validation.

## Validation status (2026-10-28)

The user-reported server results confirm:

- regenerated-input analytic covariance: **65/65 checks exactly equal**;
- catalog preprocessing, pair/model generation and signal/noise simulations:
  **146/146 checks passed, 122 exactly equal**;
- the four `offset=0` comparison panels visually reproduce the earlier figures.

The non-exact simulation differences and numerical tolerances are recorded in
the validation document. The `offset=1`, correlation-matrix and noise-increment
panels have not yet received real-data visual confirmation. Numerical
regression is not independent validation of the physical model.

## Repository layout

- `examples/xmm_nside512/`: established scientific scripts and reusable notebook.
- `docs/`: user instructions, conventions and validation record.
- `validation/`: synthetic checks and the original idark reproduction runners.
  These runners contain study-specific paths and reference results; new users
  should begin with the analytic quickstart instead.

No software license has yet been selected for this private review version.
Retain source provenance; confirm author attribution and licensing before a
public release. See [sharing notes](docs/release_and_hsc_application.md).
