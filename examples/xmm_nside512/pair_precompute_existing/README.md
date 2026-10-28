# Existing XMM Nside=512 pair-precomputation code

`pair_chunk_worker.py` and `merge_pair_chunks.py` are copied from the existing `local/XMM/Nside512/` workflow; the pair-enumeration and merge algorithms are preserved. `prep_pair_metadata.py` is also copied from that workflow. The accompanying `config_custom_10bins.py` only replaces the original machine-specific paths with environment variables and a portable output root.

The metadata stage imports the established `covariance.galaxy_catalog_to_healpix` function from the original module copied from `/lustre/work/Kousuke0128/2PCF/funcs` under `legacy_2pcf_funcs/`. Set `SHEAR_2PCF_FUNCS` only if you intentionally want to use another copy of that same module. Set `SHEAR_PAIR_WORK_ROOT` to a fresh work directory and pass its `outputs/precompute/` subfolder explicitly with `prepare_shear.py --pair-base`. Also set:

- `SHEAR_PAIR_CATALOG`: authorized XMM FITS catalog;
- `SHEAR_PAIR_THETA_NPZ`: user-supplied bin-settings NPZ containing `theta` and `theta_edges` in degrees;
- `SHEAR_PAIR_WORK_ROOT`: optional output directory (defaults beside these scripts).

Then run `python prep_pair_metadata.py`, submit the pair workers with `bash submit_pair_chunks.sh`, and run `python merge_pair_chunks.py` after every chunk succeeds. The submission wrapper preserves the original 2,000-active-pixel chunking and the PBS defaults (tiny queue, 1 CPU, 8 GB, 48 hours); set `SHEAR_PAIR_CHUNK_SIZE`, `SHEAR_PAIR_QUEUE`, `SHEAR_PAIR_MEMORY`, or `SHEAR_PAIR_WALLTIME` only when adapting the resource request. It uses `SHEAR_PAIR_WORK_ROOT/outputs/precompute/` for the metadata and generated pair files, matching the path expected by `prepare_shear.py`.

The metadata stage preserves the existing weighted HEALPix mapping via `galaxy_catalog_to_healpix`. It writes the pair metadata needed by the existing pair-worker. It does not calculate the unused scalar `xi_table.npz` side product. Shear `xi_+` and `xi_-` are computed separately by `shear_core.py:xi_from_cl` with spin-2 functions. No pair-generation or correlation algorithm has been replaced.

For a catalog-to-result run, generate pair products with these scripts. Historical products are optional regression references, not required inputs. The serial worker/merge sequence is included in `docs/analytic_quickstart.md`; PBS submission is optional.
