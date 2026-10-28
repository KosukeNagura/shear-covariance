# Existing 2PCF utility used for pixel metadata

`covariance.py` is copied from the user's existing idark installation at
`/lustre/work/Kousuke0128/2PCF/funcs/covariance.py`. The pair-metadata stage
uses its established `galaxy_catalog_to_healpix` function to create the
weighted HEALPix map. No replacement implementation is supplied. During the sharing review, unused
commented example calls were removed; all executable code and function
docstrings remain unchanged.

The scalar `w_theta` utility is intentionally not included: it only produced a
legacy `xi_table.npz` side product, which is not consumed by shear preparation,
the shear correlation calculation, or covariance calculation. Shear
`xi_+` and `xi_-` are calculated by `shear_core.py:xi_from_cl` with spin-2
spherical functions.

Before publishing this repository, confirm the source file's redistribution
terms and add any required attribution or license information.
