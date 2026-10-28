# Private review and future sharing

## Selected destination

The user selected `KosukeNagura/shear-covariance`, **Private**, on 2026-10-28.
The repository has been created. This source snapshot is prepared for its
first upload; publication/upload completion is recorded separately after the
GitHub file tree is verified. No collaborators have been invited and no
software license has been selected.

## Included material

The repository contains runnable scientific scripts, developer tests,
documentation, two comparison notebooks and study-specific validation runners.
Survey catalogs, spectrum archives, pair products, simulation outputs, logs,
checkpoints and saved notebook outputs are excluded. The shared source keeps
the original directory structure so sibling-module imports continue to work.

The reusable notebook is under `examples/`. The notebook under `validation/`
has the original server paths filled in for the reproduction record; other
users should use the reusable notebook and select their own output files.
The validation PBS scripts also retain idark paths for that record. They are
not portable launch commands for another compute cluster.

## Scientific scope

The tested setup is XMM, one HSC Y3 region, tomographic bin 4, Nside512,
ten angular bins and an EE-only spectrum. The documented workflow computes
thinned signal-only and noise-added covariance. Existing full-mode functions
remain in the source but unthinned production is not required or part of this
release validation. No B-mode null-test application is claimed.

Analytic covariance uses grid pixel pairs; comparison simulations use catalog
TreeCorr galaxy pairs. Both numerical reproduction jobs passed. All four
`offset=0` panels were visually checked against the historical figures.
Real-data visual review of the remaining panels is pending. See
[the acceptance record](reproducibility_validation.md) for the detailed scope.

## Editorial changes after the successful server runs

README and user instructions were reorganized around declared inputs and
outputs. Essential units, data-vector ordering, per-component noise variance,
HT normalization and the mixed-plus-noise increment are explained. Obsolete
commented example calls for unrelated scalar routines were removed from the
inherited utility; all functions remain intact. Comments clarify same-pixel
limits and matrix indexing. No executable Python statement or docstring was
changed: AST equality was checked for all 22 Python files before and after
this editorial pass. Eight small validation tests passed afterward.

## Before public distribution or a new survey application

Confirm source authorship/attribution and choose the software license with the
user and relevant authors. A private review upload does not settle those
questions. Catalog access and any data distribution remain separate from the
code repository. For a new sample, explicitly adapt and validate the catalog
schema/selection, weights, angular bins, signal and noise conventions. The
current guards intentionally prevent unsupported configurations from being
silently treated as validated.
