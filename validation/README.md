# Catalog-to-simulation reproduction (2026-10-27)

This validation starts independently of covariance-only job **605073**. It
regenerates pixel metadata, all ten pair lists, the shear model, and signal and
shape-noise simulations with the packaged scientific scripts. It compares each
stage with the historical result and stops at the first discrepancy.

This is not yet the final figure-reproduction result. After this job passes,
the new `model.npz` must feed the thinned analytic calculation; shape noise
can be calibrated directly from the catalog. The notebook then reads the new
analytic output and `noise_3000sims.npz`. Neither
unthinned jobs nor changes to the existing production environment are needed.

## Input and reference roles

Generation inputs:

- The original authorized XMM FITS catalog.
- The original NPZ spectrum, key `cl_theory` (a theoretical spectrum saved in
  the convergence simulation archive, not the measured simulated spectrum).
- The original angular-bin definitions and interpolation-grid size. The runner
  exports only `theta` and `theta_edges` to `bin_settings.npz` and records the
  size of `theta_grid`; it does not copy model values into the generated model.
- Nside512, bin4, seeds starting at 0, and the existing scientific conventions.

Comparison references, never passed to the generation subprocesses:

- Original model and its pair-file paths/checksums, for numeric comparison.
- Original signal and noise simulation results, for matched-seed comparison.

The original metadata/pair workers, merger, model builder and signal simulator
are called directly. No replacement pair-construction or covariance algorithm
has been introduced. The only simulation change is making historical
reference checks optional. The catalog noise calibration, random generation,
map generation and TreeCorr measurement are preserved. The synthetic integration
test checks that omitting references leaves every scientific output identical.

## Execution order

1. Record code/input hashes and dependency versions; verify catalog/spectrum identity.
2. Run the small synthetic simulation regression, including with/without-reference equality.
3. Rebuild metadata from FITS; compare pixel IDs, weights and bin settings.
4. Run every original pair chunk (size 2000), check all chunk files exist, merge,
   and require exact pair-array ordering/content and counts.
5. Build the model from the regenerated pair lists and original input spectrum;
   compare numeric arrays including the signal covariance and predicted xi.
6. Generate 3 signal and noise realizations independently; compare against the
   corresponding historical seeds (including covariance of that subset).
7. Only if all earlier checks pass, generate 3000 signal and noise realizations
   and compare per-realization measurements and sample covariances.

Exact equality is recorded separately from tolerance-based acceptance. Pair
indices and integer metadata must match exactly. Model floating comparisons
use `rtol=1e-12, atol=1e-18`; simulation xi uses the existing signal-regression
tolerance `rtol=1e-9, atol=1e-14`. Simulation covariance comparison uses
`rtol=2e-9, atol=1e-20`, not NumPy's much larger default absolute tolerance.
These are roundoff allowances, not sample-variance allowances: seeds match.
If a tolerance is exceeded, inspect the recorded difference before changing it.

## Paths

Local package:
`/Users/kosuke/.codex/.chatgpt-projects/g-p-691fd977acb88191a570d4fb2f6a97b4/shear-covariance-reproducibility`

Uploaded package:
`/lustre/work/Kousuke0128/local/lensing/shear-covariance-reproducibility_validation`

New validation products (separate from job 605073):
`/lustre/work/Kousuke0128/local/lensing/shear-covariance-reproducibility_validation_from_catalog`

The PBS example is submitted from the uploaded package root. It creates a
Python-source-only snapshot `code.JOB_ID` and a fresh `run.JOB_ID` directory.
It neither writes to historical results nor overwrites job 605073's files.
The `log/from_catalog.JOB_ID.txt` file receives live stdout/stderr directly on
the shared filesystem, including child-process output and `/usr/bin/time`.

`run.JOB_ID/validation_report.json` records commands, exit codes, input/code
hashes, versions and comparisons. Its final success status is
`upstream_passed`, explicitly not full-pipeline/figure reproduction. Partial
failure reports are retained. Re-running the runner against an existing work
directory is refused to prevent stale products from silently passing.

## Submission on idark (user operated)

After uploading the updated simulation files and `validation/`, activate the
existing `cellcov_3` environment, change to the uploaded package root, and run:

```sh
sha256sum -c --quiet validation/source_sha256.txt
python validation/test_validation_checks.py
bash -n validation/run_from_catalog.pbs
qstat -u "$USER"
qstat -Qf mini
qsub -q mini validation/run_from_catalog.pbs
```

The example requests 8 CPUs, 32 GB and a 72-hour limit. These are resource
requests, not a completion-time forecast; the regenerated preprocessing time
has not yet been measured. Do not delete other jobs merely because this one is
queued. Inspect this job's scheduler comment if it waits.

## Local checks completed and remaining limits

- Python syntax and PBS shell syntax checked.
- Five direct comparison/noise-calibration checks pass with NumPy and stdlib.
- Four existing NumPy/SciPy mathematical tests pass by direct invocation (Wick/quadratic trace, exact HT expectation, Jacobi recurrence/zero lag, and four shear-block expressions); the pytest-dependent spherical-harmonic test was not run locally.
- Original noise calibration, catalog loader and noise RNG function ASTs match
  the original local scripts; `shear_core.py` is byte-identical to the original.
- Original model builder and signal simulator differences were inspected:
  they concern input-path configuration only.
- Local Python lacks Healpy, Astropy and TreeCorr. Consequently the new full
  synthetic simulation integration test and real-catalog regeneration are
  **not yet executed here**; they run on the server in the existing environment.
- No server connection or job submission is performed by the assistant.
