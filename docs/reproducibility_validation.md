# Reproduction validation plan

Run this only after the package contents and configuration are frozen. Use the
same authorized XMM HSC Y3 catalog, bin-4 and B-mode-mask selection, spectrum,
theta bins and random seeds. Regenerate pair products, model and simulations; compare them against historical products as the existing analysis. Record
the package revision, Python and dependency versions, input paths and SHA-256
hashes, and each command's exit status. Never overwrite the reference outputs.

## Stage gates

1. **Environment and static checks.** Record the existing environment without
   installing or updating its dependencies. Use a separate environment only when
   deliberately testing a clean installation. Run the unit tests and syntax checks. Confirm the
   scripts refer only to declared inputs and that the package has no embedded
   catalog, credentials, or large generated outputs.
2. **Catalog selection and pixel metadata.** Run the included metadata script
   on the authorized catalog. Check the selected object count against the
   recorded XMM bin-4 reference, then compare `pair_meta.npz` values (`nside`,
   active pixel IDs, weights, bin edges) against the previous metadata. Pixel
   IDs and bin definitions should match exactly; floating weights must match
   within a documented numerical tolerance.
3. **Pair-list reproduction.** Run the existing pair worker and merger. Check
   per-bin pair counts and compare each merged pair array to the saved original
   arrays. Counts and integer pair indices should match exactly. Do not
   continue to covariance evaluation if they differ; first identify whether
   the input selection, HEALPix indexing, ordering, or script version changed.
4. **Shear model reproduction.** Run `prepare_shear.py` from the reproduced
   pair lists with the original spectrum. Compare model inputs, hashes,
   counts, and xi arrays with the saved `model.npz`; require exact agreement
   where the same deterministic operations are used, otherwise document and
   justify a numerical tolerance before comparing.
5. **Catalog simulation reproduction.** Run a small fixed-seed pilot, then the
   configured realization count. Check that the signal-only part agrees with
   the stored signal reference for each matched seed, and that signal, noise,
   and total outputs are finite and have the expected shapes. Compare the
   final catalog covariance with the prior result using the identical
   TreeCorr settings and realization seeds.
6. **Analytic covariance reproduction.** Run the thinned calculation with its
   original keep probability, realization count, and seed. Check completion
   metadata and finite arrays, then compare it with the original completed
   thinned result. Include the estimator's sampling uncertainty. The
   non-thinned calculation is outside the scope of this validation.
7. **Notebook reproduction.** Restart the notebook kernel and run all cells
   from top to bottom using the regenerated simulation and covariance files.
   Confirm plotted data ranges, angle bins, units, covariance block ordering,
   and error bars match the documented figure definitions.

## Acceptance record

For each stage, save: command, input hashes, output hashes, software versions,
exit status, and a short pass/fail note. Record exact checks separately from
floating-point comparisons and state every tolerance explicitly. If the
catalog TreeCorr measurement and grid covariance differ, report that as an
estimator comparison; it does not by itself mean either implementation failed
to reproduce its own reference result.

Do not make a public reproducibility claim until all applicable stage gates
pass. If a required historical input or reference output is unavailable, mark
that gate as not verified rather than substituting a different input.

## Current server validation scope (2026-10-27)

Job 605073 reruns only the thinned covariance using historical model and simulation inputs.
It does not validate catalog preprocessing, model construction or simulation generation.
`validation/run_from_catalog.py` and its PBS example independently regenerate and compare
those upstream stages, including a 3-realization pilot followed by 3000 realizations.
They need not wait for job 605073. See `validation/README.md` for exact scope and paths.
Even after both jobs pass, run the thinned covariance on the regenerated upstream outputs
and execute the notebook using those regenerated products before claiming full plot reproduction.

## Covariance-only result reported by the user (2026-10-27)

Job 605073 completed 100/100 thinnings in 38:31.60 with exit status 0.
The user's server-side comparison reports exact array equality for all 20
listed settings/result fields, including covariance_joint, covariance_sem,
covariance_blocks, covariance_realizations, bin_seeds, retained_counts and
pixel_noise_variance. This establishes reproduction of the thinned covariance
with the historical model/catalog/simulation inputs. The three stored
simulation_covariance arrays are copied inputs, not independently regenerated
simulation results in this job.

The archive is not identical in all fields: covariance_noise_added is (20,20)
in the historical output and empty (0,) in the package output. Code inspection
shows this auxiliary array is total HT covariance minus a signal-only HT result
at a package-relative path, and is empty when that result is absent. The notebook
therefore skips its final auxiliary noise-increment panels; total covariance
comparison plots use covariance_joint/covariance_sem and do not require it.
This auxiliary-output difference was subsequently resolved by job 605209,
using an explicit matching regenerated signal-only HT result (see below). It does not indicate
missing shape noise in the total covariance, which matched exactly.

The user subsequently reported successful completion of job 605199; the
2026-10-28 acceptance record below supersedes its earlier running status.
No server access was made by the assistant.


## Analytic-first follow-up

The user prioritizes analytic verification over completion of the long
simulation runs. `validation/run_analytic_validation.py` verifies the already
regenerated metadata/pair/model evidence and computes both HT covariances from
that new model, using catalog-derived noise normalization directly. It has no
simulation-result dependency and need not wait for job 605199. Its explicit
signal-only reference also fills and verifies the noise-increment auxiliary
array. Job 605209 has now passed all 65 comparisons exactly, including this
auxiliary array. See `validation/ANALYTIC.md` for execution scope and the record
below for completed results. Real-data notebook execution remains pending.

## Review snapshot preparation

The comparison notebook now accepts explicit absolute paths for both regenerated
outputs and prints the paths it reads. Its calculation and plotting cells retain
the existing covariance definitions, arcmin conversion and plotting style. The
3000-realization requirement prevents a pilot being mislabeled as a production
comparison. This input-selection change is local preparation; it does not change
the code snapshots used by the running server jobs.

A source archive and its file manifest can be prepared before server completion.
Archive integrity, local tests and synthetic notebook execution do not establish
real-catalog figure reproduction. The final server analytic report and eventual
regenerated simulation comparisons must be recorded separately.

Local checks on this review snapshot: the five validation-comparison tests and
three analytic-input tests passed. The notebook's code cells executed and four
figures rendered with synthetic input arrays; mismatched model hashes and a
3-realization pilot were rejected. No new dependency was installed. This checks
execution and input handling only, not agreement of real-catalog figures.

## Server acceptance results reported 2026-10-28

Evidence: the user supplied both JSON report summaries, the 24 non-exact
comparison records, and both final job logs. No direct server access was used.

| Job | Scope | Status | Comparisons | Wall time / exit |
| --- | --- | --- | --- | --- |
| 605209 | New-model signal-only and noisy HT, plus paired noise increment; no simulation input | analytic_passed | 65/65 exact | 1:41:42 / 0 |
| 605199 | Catalog metadata, full pair arrays, model, pilot and production signal/noise simulations | upstream_passed | 146/146 passed, 122 exact | 7:38:16 / 0 |

The 24 non-exact checks are simulation-derived correlation/covariance arrays:
12 pilot and 12 production checks. Across all 3000-realization production
checks, the largest absolute difference in correlations (including imaginary
components) is 5.827586677109586e-18; in covariances it is
3.7223127564886245e-24. Their absolute tolerances were 1e-14 and 1e-20,
respectively (relative tolerances 1e-9 and 2e-9). Pilot maxima were
2.7511630126819675e-18 and 5.469731800506896e-23. All pass even using the
absolute tolerance alone. The reports do not identify the mechanism causing
these tiny differences; no causal claim is made here.

This establishes numerical regression reproduction for the tested XMM setup,
not independent proof of the physical model or general applicability to other
catalogs. Both standalone job logs describe work outside their own scope as
pending; the two reports together confirm the analytic and simulation stages.
The remaining figure check uses
`validation/notebooks/compare_regenerated_shape_noise.ipynb`, which has both
regenerated input paths set. Plotting and covariance/error definitions are
unchanged from the existing notebook. At that stage real-data figure execution was pending. The subsequent
offset=0 visual review is recorded below.

## Visual review reported 2026-10-28

The user supplied new and historical screenshots of all four offset=0 panels.
Their points, curves and error bars show no visible change; the earlier small-
angle grid/catalog differences and the cross-block endpoint behavior are
reproduced. This supports the numerical regression results. This is a visual
check of offset=0 only; offset=1, correlation-matrix and noise-increment panels
have not yet been visually confirmed with the real outputs.
