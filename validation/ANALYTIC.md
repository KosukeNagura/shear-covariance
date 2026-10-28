# Prioritize analytic reproduction independently of simulations

The user's priority is to validate the primary analytic calculation first and
prepare GitHub sharing if that succeeds. The 3000-realization simulation check
can continue separately and is no longer a gate for this analytic test.

## What the new job does

`run_analytic_validation.py` reads the completed metadata, pair and model check
records from `run.605199.idark/validation_report.json`. It does not require that
job to finish, and does not read either pilot or production simulation arrays.
It verifies the regenerated model and pair files again, calibrates noise directly
from the catalog, then computes signal-only and signal-plus-noise HT covariances
with p=0.001, 100 realizations and base seed 20260828. Both use the regenerated
`model.npz` and pair files; historical HT outputs are used only for comparison.

The noise calculation now accepts `--signal-only-reference` to save its difference
from the newly computed signal-only HT output. It checks matching model hash,
seeds and retained counts. The difference is compared with the historical
`covariance_noise_added`, resolving the missing auxiliary array from job 605073.
This difference includes mixed and noise-only contributions; it is not the
catalog simulation's noise-only covariance.

No simulation estimate is used to derive the noise amplitude: the established
`load_inputs` / `noise_parameters` code computes it from catalog moments. The
new analytic CLI defaults to this direct route. Supplying `--sims` explicitly
retains the former archive-input route. No Wick kernel, pair thinning, geometry,
normalization or RNG implementation has been changed.

In direct mode the output records `noise_calibration_source=catalog`,
`simulation_nsims=0` and empty `simulation_covariance_*` arrays. These are absent
comparison data, not estimated zero matrices. The separate comparison notebook
reads its simulation covariance from an independently supplied simulation file.

## User-operated submission

Upload the updated package to the same server package directory. Job 605199
already runs from its own code snapshot, so this upload does not change its code.
From the uploaded package root in the existing `cellcov_3` environment:

```sh
sha256sum -c --quiet validation/source_sha256.txt &&
python validation/test_analytic_inputs.py &&
bash -n validation/run_analytic_validation.pbs &&
python validation/run_analytic_validation.py --check-ready \
  --upstream-run /lustre/work/Kousuke0128/local/lensing/shear-covariance-reproducibility_validation_from_catalog/run.605199.idark &&
qsub -q mini2 validation/run_analytic_validation.pbs
```

Do not add an `afterok:605199` dependency: the prerequisites needed here have
already been created and checked, and the readiness check verifies that fact.
The job requests 1 CPU, 16 GB and 24 hours, runs the two HT calculations serially,
and uses fresh per-job output and code snapshot directories under
`/lustre/work/Kousuke0128/local/lensing/shear-covariance-reproducibility_validation_analytic`.
Live log: `log/analytic.JOB_ID.txt`. Report: `run.JOB_ID/analytic_report.json`.
The success status is `analytic_passed`; full simulation/figure validation stays
explicitly pending. Both numeric agreement and exact equality are reported.

## Checks already run locally

The three direct tests in `test_analytic_inputs.py` pass. They check the matched
subtraction guards, rejection of incomplete upstream evidence, and equality of
the actual covariance kernels' outputs through the two calibration-input routes.
Only FITS/HEALPix I/O is mocked in the latter small four-pixel numerical fixture;
it does not establish real-catalog reproduction. The real-catalog job is the
required integration check. The five existing NumPy/stdlib checks also pass.
No new software installation or server connection was made by the assistant.
