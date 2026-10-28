import os
import sys
import numpy as np
import healpy as hp
from astropy.io import fits

DEFAULT_FUNCS_PATH = os.path.join(os.path.dirname(__file__), "legacy_2pcf_funcs")
FUNCS_PATH = os.environ.get("SHEAR_2PCF_FUNCS", DEFAULT_FUNCS_PATH)
if FUNCS_PATH not in sys.path:
    sys.path.insert(0, FUNCS_PATH)
from covariance import galaxy_catalog_to_healpix

from config_custom_10bins import (
    IZ, NSIDE, AREA, WORK_ROOT, CATALOG_PATH,
    THETA_DEG, THETA_EDGES_DEG, THETA_EDGES_RAD, NBINS,
    ensure_dirs, pair_meta_path
)

def main():
    ensure_dirs()

    print("=== prep_pair_metadata.py started ===", flush=True)
    print(f"WORK_ROOT = {WORK_ROOT}", flush=True)
    print(f"AREA      = {AREA}", flush=True)
    print(f"NSIDE     = {NSIDE}", flush=True)
    print(f"NBINS     = {NBINS}", flush=True)
    print(f"THETA_EDGES_DEG = {THETA_EDGES_DEG}", flush=True)

    # -------------------------
    # read XMM catalog
    # -------------------------
    with fits.open(CATALOG_PATH) as hdul:
        data = hdul[1].data

    mask = (data["b_mode_mask"] == True) & (data["hsc_y3_zbin"] == IZ)

    ra_cat = np.asarray(data["i_ra"][mask], dtype=float)
    dec_cat = np.asarray(data["i_dec"][mask], dtype=float)
    w_cat = np.asarray(data["i_hsmshaperegauss_derived_weight"][mask], dtype=float)

    print(f"selected galaxies = {len(ra_cat)}", flush=True)

    # -------------------------
    # make HEALPix map
    # -------------------------
    count_map, weight_map, delta_map = galaxy_catalog_to_healpix(
        ra_deg=ra_cat,
        dec_deg=dec_cat,
        nside=NSIDE,
        weight=w_cat
    )

    print("weight_map constructed.", flush=True)
    print(f"nonzero pixels = {np.count_nonzero(weight_map)}", flush=True)

    # -------------------------
    # metadata for pair search
    # -------------------------
    pix_idx = np.where(weight_map != 0.0)[0].astype(np.int32)
    weights = weight_map[pix_idx].astype(np.float32)

    theta_pix, phi_pix = hp.pix2ang(NSIDE, pix_idx)
    theta_vecs = np.array(hp.ang2vec(theta_pix, phi_pix), dtype=np.float32)

    npix = hp.nside2npix(NSIDE)
    pix_to_local = np.full(npix, -1, dtype=np.int32)
    pix_to_local[pix_idx] = np.arange(len(pix_idx), dtype=np.int32)

    np.savez(
        pair_meta_path(),
        iz=IZ,
        nside=NSIDE,
        area=AREA,
        theta=THETA_DEG,
        theta_edges=THETA_EDGES_DEG,
        theta_edges_rad=THETA_EDGES_RAD,
        pix_idx=pix_idx,
        weights=weights,
        theta_vecs=theta_vecs,
        pix_to_local=pix_to_local,
    )
    print(f"saved pair meta to {pair_meta_path()}", flush=True)

    print("=== prep_pair_metadata.py finished ===", flush=True)

if __name__ == "__main__":
    main()
