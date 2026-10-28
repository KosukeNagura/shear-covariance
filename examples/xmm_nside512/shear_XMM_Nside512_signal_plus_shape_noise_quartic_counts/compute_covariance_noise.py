#!/usr/bin/env python3
"""Shear Gaussian covariance with catalog-calibrated pixel shape noise.

Uses the existing signal covariance matrix S and adds the variance of the
weighted per-pixel mean of independent catalog galaxy noise on S[u,u]. The
same quartic-counts kernel and exact HT self-pair correction as signal-only
are retained. Output is the grid estimator's closed total covariance.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SIGNAL = HERE.parent / 'shear_XMM_Nside512_signal_only_quartic_counts'
sys.path.insert(0, str(SIGNAL))

from shear_core import (file_hash, exclusive_save, checkpoint_save, pair_geometry,
                        bernoulli_pairs, identical_sum, raw_pair_sum, assemble_blocks)

DEFAULT_MODEL = SIGNAL / 'outputs/model.npz'
DEFAULT_CATALOG = HERE.parent / 'data' / 'XMM.fits'


def read_pixel_noise(catalog, model, sigma_component_squared, expected_ngal):
    import healpy as hp
    from astropy.io import fits

    nside, zbin = int(model['nside']), int(model['zbin'])
    with fits.open(catalog, memmap=True) as hdul:
        table = hdul[1].data
        sel = (table['b_mode_mask'] == True) & (table['hsc_y3_zbin'] == zbin)
        ra = np.asarray(table['i_ra'][sel], dtype=np.float64)
        dec = np.asarray(table['i_dec'][sel], dtype=np.float64)
        w = np.asarray(table['i_hsmshaperegauss_derived_weight'][sel], dtype=np.float64)
    good = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(w) & (w > 0)
    ra, dec, w = ra[good], dec[good], w[good]
    if len(w) != expected_ngal:
        raise ValueError(f'Catalog selection has {len(w)} galaxies; simulation has {expected_ngal}')
    pix = hp.ang2pix(nside, ra, dec, lonlat=True)
    npix = hp.nside2npix(nside)
    sumw = np.bincount(pix, weights=w, minlength=npix)[model['pix_idx']]
    sumw2 = np.bincount(pix, weights=w*w, minlength=npix)[model['pix_idx']]
    if not np.allclose(sumw, model['weights'], rtol=3e-6, atol=1e-5):
        raise ValueError('Catalog pixel weights do not match signal-only model')
    if np.any(sumw <= 0):
        raise ValueError('Active model pixel has no catalog weight')
    # Variance of weighted mean: Var(sum w_g n_g / sum w_g).
    return sigma_component_squared * sumw2 / sumw**2


def catalog_noise_parameters(catalog, model):
    """Use the established simulation loader, without generating a simulation."""
    noise_dir = HERE.parent / 'shear_XMM_Nside512_signal_plus_shape_noise_sims'
    sys.path.insert(0, str(noise_dir))
    from simulate_shear_shape_noise import load_inputs
    ra, dec, w, pix, parameters, names = load_inputs(
        catalog, int(model['nside']), int(model['zbin']))
    return parameters, len(w)


def checked_noise_increment(path, joint, seeds, counts, probability, base_seed, model_hash=None):
    """Subtract only a signal-only result with the same thinning realization."""
    with np.load(path) as f:
        base_joint = np.asarray(f['covariance_joint'], dtype=np.float64)
        same = (float(f['keep_probability']) == probability and int(f['base_seed']) == base_seed
                and np.array_equal(f['bin_seeds'], seeds)
                and np.array_equal(f['retained_counts'], counts))
        if model_hash is not None:
            same = same and str(f['model_sha256']) == model_hash
    if not same:
        raise ValueError('Signal-only result has different model or thinning realization; refusing subtraction')
    if base_joint.shape != joint.shape or not np.all(np.isfinite(base_joint[np.isfinite(joint)])):
        raise ValueError('Signal-only result has incompatible shape or incomplete calculated cells')
    # This increment contains both signal-noise mixed and pure-noise terms.
    return joint - base_joint


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--model', type=Path, default=DEFAULT_MODEL)
    ap.add_argument('--sims', type=Path, help='Optional legacy calibration/simulation input; omit to calibrate directly from the catalog')
    ap.add_argument('--signal-only-reference', type=Path, help='Matching signal-only HT output used only for the saved noise increment')
    ap.add_argument('--catalog', type=Path, default=DEFAULT_CATALOG)
    ap.add_argument('--mode', choices=['full', 'thinned'], required=True)
    ap.add_argument('--row', type=int, help='full mode: calculate only this upper-triangle row')
    ap.add_argument('--realizations', type=int, default=100)
    ap.add_argument('--base-seed', type=int, default=20260828)
    ap.add_argument('--left-chunk', type=int, default=0)
    ap.add_argument('--right-chunk', type=int, default=2000)
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    if args.mode == 'thinned' and args.row is not None:
        raise ValueError('--row only applies to full mode')
    if args.realizations < 1 or args.left_chunk < 0 or args.right_chunk < 1:
        raise ValueError('Invalid realization or chunk size')
    if args.mode == 'full' and args.row is not None and not 0 <= args.row < 10:
        raise ValueError('--row must be in 0..9')

    with np.load(args.model) as f:
        model = {k: f[k] for k in f.files}
    if args.signal_only_reference is not None and not args.signal_only_reference.is_file():
        raise FileNotFoundError(args.signal_only_reference)
    if args.sims is None:
        # Reuse the original catalog selection and noise normalization verbatim.
        noise_meta, ngal = catalog_noise_parameters(args.catalog, model)
        sigma = noise_meta['sigma_component_squared']
        sim_total = sim_signal = sim_noise = np.array([])
        nsims = 0
    else:
        with np.load(args.sims) as f:
            sigma = float(f['sigma_component_squared'])
            ngal = int(f['n_galaxies'])
            theta = np.asarray(f['theta'])
            sim_total = np.asarray(f['covariance_total'])
            sim_signal = np.asarray(f['covariance_signal'])
            sim_noise = np.asarray(f['covariance_noise'])
            nsims = int(f['nsims']) if 'nsims' in f else -1
            nside, zbin = int(f['nside']), int(f['zbin'])
        if (nside != int(model['nside']) or zbin != int(model['zbin']) or
                not np.allclose(theta, model['theta'], rtol=0.05)):
            raise ValueError('Simulation and grid model geometry differ')

    variance = read_pixel_noise(args.catalog, model, sigma, ngal)
    # Independent isotropic shear noise has covariance v_i I at one pixel and
    # zero covariance between distinct pixels, in the fixed local basis.
    total_signal = np.asarray(model['signal'], dtype=np.float64).copy()
    pix = np.arange(len(variance))
    total_signal[pix, pix, 0, 0] += variance
    total_signal[pix, pix, 1, 1] += variance
    print(f'galaxies={ngal} active_pixels={len(variance)} sigma_component^2={sigma:.10g} '
          f'mean_grid_noise_variance={np.mean(variance):.10g}', flush=True)

    pairs = []
    for b, (path, digest) in enumerate(zip(model['pair_paths'], model['pair_hashes'])):
        pth = Path(path)
        if file_hash(pth) != digest:
            raise ValueError(f'Pair list changed in bin {b}: {pth}')
        pairs.append(np.load(pth, mmap_mode='r'))

    probability = 1.0 if args.mode == 'full' else 0.001
    nrep = 1 if probability == 1 else args.realizations
    default_name = ('full' if args.row is None else f'full_row_{args.row:02d}') if probability == 1 else 'ht_p0.001'
    outroot = HERE if probability == 1 else HERE.parent / 'shear_XMM_Nside512_signal_plus_shape_noise_quartic_counts_down_sampling'
    output = args.output or outroot / 'outputs' / f'{default_name}.npz'
    if output.exists():
        raise FileExistsError(f'Refusing to overwrite existing output: {output}')
    model_hash = file_hash(args.model)
    config = json.dumps(dict(model_sha256=model_hash, mode=args.mode, row=args.row,
                             probability=probability, nrep=nrep, seed=args.base_seed,
                             catalog=str(args.catalog), sims=str(args.sims)), sort_keys=True)
    checkpoint = output.with_suffix('.checkpoint.npz')
    blocks = np.full((nrep, 10, 10, 2, 2), np.nan)
    counts = np.zeros((nrep, 10), dtype=np.int64)
    seeds = np.zeros((nrep, 10), dtype=np.uint64)
    completed = 0
    if checkpoint.exists():
        with np.load(checkpoint) as f:
            if str(f['config']) != config:
                raise ValueError('Checkpoint configuration differs; use a new output path')
            completed = int(f['completed']); blocks[:] = f['blocks']
            counts[:] = f['retained_counts']; seeds[:] = f['bin_seeds']
        print(f'Resuming realization {completed + 1}/{nrep}', flush=True)

    started = time.monotonic()
    seqs = np.random.SeedSequence(args.base_seed).spawn(nrep)
    den = np.outer(model['full_weight_sums'], model['full_weight_sums'])
    for r in range(completed, nrep):
        selected = []
        for b, seq in enumerate(seqs[r].spawn(10)):
            seeds[r, b] = seq.generate_state(1, dtype=np.uint64)[0]
            selected.append(bernoulli_pairs(pairs[b], probability, np.random.default_rng(seq)))
            counts[r, b] = len(selected[-1])
        geom = [pair_geometry(x, model['vectors'], model['weights']) for x in selected]
        diag = [identical_sum(g, total_signal) for g in geom] if probability < 1 else None
        rows = range(10) if args.row is None else [args.row]
        for a in rows:
            for b in range(a, 10):
                if np.all(np.isfinite(blocks[r, a, b])):
                    continue
                raw = raw_pair_sum(geom[a], geom[b], total_signal,
                                   args.right_chunk, args.left_chunk)
                # Distinct retained pairs have inclusion probability p^2.
                # The identical pair has probability p, requiring the diagonal correction.
                if probability < 1:
                    raw = (raw - (1-probability)*diag[a] if a == b else raw) / probability**2
                block = raw / den[a, b]
                blocks[r, a, b] = block
                blocks[r, b, a] = block.T
                checkpoint_save(checkpoint, config=config, completed=r,
                                blocks=blocks, retained_counts=counts, bin_seeds=seeds)
                print(f'realization={r+1}/{nrep} bin=({a},{b}) '
                      f'elapsed={time.monotonic()-started:.1f}s', flush=True)
        checkpoint_save(checkpoint, config=config, completed=r+1,
                        blocks=blocks, retained_counts=counts, bin_seeds=seeds)

    mean = np.mean(blocks, axis=0)
    sem = np.std(blocks, axis=0, ddof=1)/np.sqrt(nrep) if nrep > 1 else np.zeros_like(mean)
    joint = assemble_blocks(mean)
    if args.row is None and (not np.all(np.isfinite(joint)) or not np.allclose(joint, joint.T, atol=1e-20, rtol=1e-10)):
        raise AssertionError('Covariance output is incomplete or nonsymmetric')
    baseline_path = (SIGNAL / 'outputs' / f'full_row_{args.row:02d}.npz'
                     if probability == 1 and args.row is not None else
                     SIGNAL / 'outputs' / 'full.npz' if probability == 1 else
                     SIGNAL.parent / 'shear_XMM_Nside512_signal_only_quartic_counts_down_sampling'
                     / 'outputs' / 'ht_p0.001.npz')
    if args.signal_only_reference is not None:
        baseline_path = args.signal_only_reference
    delta = None
    if baseline_path.exists():
        delta = checked_noise_increment(baseline_path, joint, seeds, counts, probability,
                                        args.base_seed, model_hash if args.signal_only_reference is not None else None)
    output.parent.mkdir(parents=True, exist_ok=True)
    exclusive_save(output, covariance_joint=joint, covariance_sem=assemble_blocks(sem),
                   covariance_blocks=mean, covariance_realizations=np.array([assemble_blocks(x) for x in blocks]),
                   # Keep the signal-only row merger's established schema too.
                   blocks=mean, row=-1 if args.row is None else args.row,
                   covariance_noise_added=delta if delta is not None else np.array([]),
                   signal_only_reference=str(baseline_path) if delta is not None else '',
                   retained_counts=counts, bin_seeds=seeds, keep_probability=probability,
                   n_realizations=nrep, base_seed=args.base_seed,
                   model_sha256=model_hash, config=config, theta=model['theta'],
                   theta_edges=model['theta_edges'], full_weight_sums=model['full_weight_sums'],
                   full_counts=model['full_counts'], pixel_noise_variance=variance,
                   sigma_component_squared=np.float64(sigma), n_galaxies=np.int64(ngal),
                   noise_calibration_source='catalog' if args.sims is None else 'simulation_archive',
                   simulation_reference='' if args.sims is None else str(args.sims),
                   simulation_nsims=nsims,
                   simulation_covariance_total=sim_total, simulation_covariance_signal=sim_signal,
                   simulation_covariance_noise=sim_noise,
                   estimator=np.asarray('HEALPix grid pair estimator; fixed full-pair denominators'),
                   note=np.asarray('Grid-closed signal + weighted pixel shape noise; unordered i<j; no extra global factor 2'))
    print(f'Saved {output}', flush=True)


if __name__ == '__main__':
    main()
