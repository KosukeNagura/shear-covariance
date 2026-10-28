#!/usr/bin/env python3
"""Catalog-only shear signal/noise/total; Convergence Gaussian-noise convention."""
import argparse
import json
import os
import sys
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
COMMON = ROOT.parent / 'shear_XMM_Nside512_signal_only_quartic_counts'
NOISE_SEED_OFFSET = 10_000_000
CATALOG = ROOT.parent / 'data' / 'XMM.fits'
LABELS = ('signal', 'noise', 'total')


def noise_parameters(e1, e2, w):
    """Same weighted second moments/responsivities as Convergence (not centered)."""
    def mean_square(x):
        return np.sum(w * x**2, dtype=np.float64) / np.sum(w, dtype=np.float64)
    r1, r2 = 1 - mean_square(e1), 1 - mean_square(e2)
    if not (0 < r1 <= 1 and 0 < r2 <= 1):
        raise ValueError('Unphysical responsivity')
    s1, s2 = mean_square(e1 / (2*r1)), mean_square(e2 / (2*r2))
    return dict(responsivity_1=r1, responsivity_2=r2,
                sigma_g1_squared=s1, sigma_g2_squared=s2,
                sigma_component_squared=0.5*(s1+s2))


def draw_noise(ngal, variance, seed):
    # Separate RNG: never advances the global signal-map random state.
    rng = np.random.default_rng(NOISE_SEED_OFFSET + seed)
    return rng.normal(0., np.sqrt(variance), size=(2, ngal))


def load_inputs(path, nside, zbin):
    import healpy as hp
    from astropy.io import fits
    with fits.open(path, memmap=True) as h:
        t = h[1].data
        names = []
        for a in (1, 2):
            candidates = (f'i_hsmshaperegauss_e{a}', f'i_hsmshaperegauss_derived_e{a}',
                          f'i_hsmshaperegauss_shape_{a}', f'e{a}')
            names.append(next((k for k in candidates if k in t.dtype.names), None))
        if None in names:
            raise KeyError('Cannot identify ellipticity columns')
        sel = (t['b_mode_mask'] == True) & (t['hsc_y3_zbin'] == zbin)
        ra, dec, w, e1, e2 = [np.asarray(t[k][sel], float) for k in
            ('i_ra', 'i_dec', 'i_hsmshaperegauss_derived_weight', *names)]
    good = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(w) & (w > 0)
    # Do not silently drop galaxies relative to the existing signal-only catalog.
    if np.any(good & (~np.isfinite(e1) | ~np.isfinite(e2))):
        raise ValueError('Ellipticity cuts would change signal-only catalog selection')
    ra, dec, w, e1, e2 = [v[good] for v in (ra, dec, w, e1, e2)]
    pix = hp.ang2pix(nside, ra, dec, lonlat=True)
    return ra, dec, w, pix, noise_parameters(e1, e2, w), names


def main():
    import healpy as hp
    import treecorr
    sys.path.insert(0, str(COMMON))
    from shear_core import file_hash, exclusive_save, checkpoint_save
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--model', type=Path, default=COMMON/'outputs/model.npz')
    ap.add_argument('--catalog', type=Path, default=CATALOG)
    ap.add_argument('--reference-signal', type=Path, help='Optional signal regression reference; never used to generate realizations')
    ap.add_argument('--reference-convergence', type=Path, help='Optional convergence noise-calibration reference; variance is computed from the catalog')
    ap.add_argument('--nsims', type=int, default=3000)
    ap.add_argument('--seed-start', type=int, default=0)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--checkpoint-every', type=int, default=25)
    ap.add_argument('--output', type=Path)
    ap.add_argument('--preflight', action='store_true')
    ap.add_argument('--require-conda', action='store_true')
    args = ap.parse_args()
    if args.nsims < 2 or args.seed_start < 0 or args.threads < 1 or args.checkpoint_every < 1:
        raise ValueError('Invalid simulation settings')
    if args.require_conda and os.environ.get('CONDA_DEFAULT_ENV') != 'cellcov_3':
        raise ValueError('Activate cellcov_3')
    output = args.output or ROOT/f'outputs/shear_signal_shape_noise_catalog_{args.nsims}sims.npz'
    if output.exists() and not args.preflight:
        raise FileExistsError(output)
    with np.load(args.model) as f:
        m = {k: f[k] for k in ('cl_ee', 'cl_bb', 'nside', 'zbin', 'theta', 'theta_edges', 'pix_idx', 'weights')}
    nside, zbin = int(m['nside']), int(m['zbin'])
    cl = m['cl_ee']; edges = m['theta_edges']; n = len(edges)-1
    if nside != 512 or zbin != 4 or n != 10 or cl.shape != (1536,):
        raise ValueError('Expected XMM bin4 Nside512, 10 bins and ell_max=1535')
    if np.any(cl[:2] != 0) or np.any(~np.isfinite(cl)) or np.any(cl < 0) or np.any(m['cl_bb'] != 0):
        raise ValueError('Invalid EE/BB signal input')
    ra, dec, w, pix, noise_meta, e_names = load_inputs(args.catalog, nside, zbin)
    wmap = np.bincount(pix, weights=w, minlength=hp.nside2npix(nside))
    if np.count_nonzero(wmap) != len(m['pix_idx']) or not np.allclose(wmap[m['pix_idx']], m['weights'], rtol=3e-6, atol=1e-5):
        raise ValueError('Catalog/model weights differ')
    model_hash = file_hash(args.model); catalog_hash = file_hash(args.catalog)
    # Keep the successful signal-only TreeCorr settings, including default slops.
    template = treecorr.GGCorrelation(min_sep=edges[0], max_sep=edges[-1], nbins=n, sep_units='deg')
    reference = None
    if args.reference_signal is not None:
        with np.load(args.reference_signal) as f:
            reference = f['xi_catalog'].copy(); ref_seeds = f['seeds'].copy()
            ref_config = json.loads(str(f['config']))
            if str(f['model_sha256']) != model_hash or ref_config['catalog_sha256'] != catalog_hash:
                raise ValueError('Signal reference model/catalog mismatch')
            for key, value in dict(metric='Euclidean', bin_slop=template.bin_slop, angle_slop=template.angle_slop,
                                   treecorr_version=treecorr.__version__, pixwin=False, fwhm=0).items():
                if ref_config.get(key) != value:
                    raise ValueError(f'Signal reference setting differs: {key}')
        ref_rows = {int(seed): row for row, seed in enumerate(ref_seeds)}
        if not all(seed in ref_rows for seed in range(args.seed_start, args.seed_start+args.nsims)):
            raise ValueError('Signal reference does not cover requested seeds')
    if args.reference_convergence is not None:
        with np.load(args.reference_convergence) as f:
            if int(f['n_galaxies']) != len(w) or int(f['nside']) != nside or int(f['zbin']) != zbin:
                raise ValueError('Convergence reference selection differs')
            for key in ('responsivity_1', 'responsivity_2', 'sigma_g1_squared', 'sigma_g2_squared'):
                if not np.isclose(noise_meta[key], float(f[key]), rtol=1e-12, atol=0):
                    raise ValueError(f'Convergence noise normalization differs: {key}')
            if not np.isclose(noise_meta['sigma_component_squared'], float(f['sigma_kappa_noise_squared']), rtol=1e-12, atol=0):
                raise ValueError('Convergence scalar/shear component variance differs')
    config = json.dumps(dict(version=1, model_sha256=model_hash, catalog_sha256=catalog_hash,
        signal_reference_sha256=file_hash(args.reference_signal) if args.reference_signal is not None else None,
        convergence_reference_sha256=file_hash(args.reference_convergence) if args.reference_convergence is not None else None, nsims=args.nsims,
        seed_start=args.seed_start, noise_seed_offset=NOISE_SEED_OFFSET,
        noise_model='independent per-galaxy isotropic Gaussian; per-component variance matches convergence scalar',
        measurement_policy='treecorr_defaults_v2', treecorr_version=treecorr.__version__,
        numpy_version=np.__version__, healpy_version=hp.__version__, threads=args.threads,
        metric='Euclidean', bin_slop=template.bin_slop, angle_slop=template.angle_slop,
        pixwin=False, fwhm=0, **noise_meta), sort_keys=True)
    print(f'galaxies={len(w)}; config={config}', flush=True)
    print('Catalog/model correspondence: PASS; noise variance computed from catalog.', flush=True)
    print('Convergence reference: ' + ('PASS' if args.reference_convergence is not None else 'NOT REQUESTED'), flush=True)
    reference_status = 'PASS' if reference is not None else 'NOT_REQUESTED'
    print('Signal reference: ' + ('per-seed checks enabled' if reference is not None else 'NOT REQUESTED'), flush=True)
    if args.preflight:
        print('Preflight passed. No simulation or job submitted.', flush=True)
        return
    arrays = {f'xi_{label}{suffix}_sims': np.full((args.nsims, 2*n), np.nan)
              for label in LABELS for suffix in ('', '_imag')}
    cp = output.with_suffix('.checkpoint.npz'); completed = 0
    if cp.exists():
        with np.load(cp) as f:
            if str(f['config']) != config:
                raise ValueError('Checkpoint configuration mismatch')
            completed = int(f['completed'])
            if not 0 <= completed <= args.nsims:
                raise ValueError('Invalid checkpoint count')
            for key in arrays:
                arrays[key][:] = f[key]
                if not np.all(np.isfinite(arrays[key][:completed])):
                    raise ValueError('Invalid checkpoint data')
    zero = np.zeros_like(cl); started = time.monotonic()
    for row in range(completed, args.nsims):
        seed = args.seed_start+row; tick = time.monotonic()
        np.random.seed(seed)
        maps = hp.synfast([zero, cl, zero, zero], nside, lmax=len(cl)-1, new=True, pol=True, pixwin=False, fwhm=0.)
        signal = np.array([maps[1][pix], -maps[2][pix]])
        del maps
        noise = draw_noise(len(w), noise_meta['sigma_component_squared'], seed)
        fields = dict(signal=signal, noise=noise, total=signal+noise)
        for label, values in fields.items():
            gg = template.copy()
            gg.process(treecorr.Catalog(ra=ra, dec=dec, ra_units='deg', dec_units='deg', w=w,
                                       g1=values[0], g2=values[1]), num_threads=args.threads)
            arrays[f'xi_{label}_sims'][row] = np.r_[gg.xip, gg.xim]
            arrays[f'xi_{label}_imag_sims'][row] = np.r_[gg.xip_im, gg.xim_im]
        if reference is not None:
            ref = reference[ref_rows[seed]]
            if not np.allclose(arrays['xi_signal_sims'][row], ref, rtol=1e-9, atol=1e-14):
                raise ValueError(f'Signal-only regression failed at seed={seed}; no silent change permitted')
        if not all(np.all(np.isfinite(a[row])) for a in arrays.values()):
            raise ValueError(f'Nonfinite measurements at seed={seed}')
        done = row+1
        if done % args.checkpoint_every == 0 or done == args.nsims:
            checkpoint_save(cp, config=config, completed=done, **arrays)
        elapsed = time.monotonic()-started
        print(f'realization={done}/{args.nsims} seed={seed} signal_reference={reference_status} '
              f'step={time.monotonic()-tick:.2f}s elapsed={elapsed:.1f}s '
              f'ETA={elapsed/(done-completed)*(args.nsims-done):.1f}s', flush=True)
    summary = {}
    for label in LABELS:
        values = arrays[f'xi_{label}_sims']
        summary[f'xi_{label}_mean'] = values.mean(axis=0)
        summary[f'covariance_{label}'] = np.cov(values, rowvar=False, ddof=1)
    exclusive_save(output, config=config, **arrays, **summary, **noise_meta,
        theta=m['theta'], theta_edges=edges, theta_units='deg', nside=nside, zbin=zbin,
        cl_ee=cl, cl_bb=zero, nsims=args.nsims, n_galaxies=len(w),
        seeds=np.arange(args.seed_start,args.seed_start+args.nsims),
        noise_seed_offset=NOISE_SEED_OFFSET, model_sha256=model_hash, catalog_sha256=catalog_hash,
        ellipticity_columns=np.array(e_names), ordering='plus bins, then minus bins',
        estimator='catalog TreeCorr defaults; nearest-pixel Q,-U without component transport')
    print(f'Saved {output}', flush=True)


if __name__ == '__main__':
    main()
