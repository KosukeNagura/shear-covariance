#!/usr/bin/env python3
"""Prepare immutable model using the existing Nside512 convergence pair lists."""
import argparse
from pathlib import Path
import numpy as np
import healpy as hp
from shear_core import exclusive_save, file_hash, xi_from_cl, signal_covariance, pair_geometry

ROOT = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pair-base', type=Path, default=ROOT/'outputs/precompute')
    ap.add_argument('--cl-file', type=Path, required=True,
                    help='NPZ spectrum input; its key is selected with --cl-key')
    ap.add_argument('--cl-key', default='cl_theory')
    ap.add_argument('--output', type=Path, default=ROOT/'outputs/model.npz')
    ap.add_argument('--theta-samples', type=int, default=60001)
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    if args.theta_samples < 30001:
        raise ValueError('Use at least 30001 angle samples for production')
    meta_path = args.pair_base/'pair_meta.npz'
    with np.load(meta_path) as m:
        pix = m['pix_idx'].astype(int)
        weights = m['weights'].astype(float)
        nside, zbin = int(m['nside']), int(m['iz'])
        edges = m['theta_edges'].astype(float)
    if nside != 512 or zbin != 4 or len(edges) != 11:
        raise ValueError('Expected Nside512, zbin4, ten bins')
    if not np.all(np.isfinite(weights) & (weights > 0)):
        raise ValueError('Invalid pixel weights')
    # Reconstruct exact float64 centres: old metadata stores float32 vectors.
    vectors = np.array(hp.pix2vec(nside, pix)).T
    lmax = 3*nside - 1
    with np.load(args.cl_file) as f:
        spectrum = np.asarray(f[args.cl_key], float)
    if spectrum.ndim != 1 or len(spectrum) <= lmax:
        raise ValueError('Input spectrum must be C_ell indexed by ell=0..lmax')
    cl = spectrum[:lmax+1].copy()
    cl[:2] = 0
    if np.any(~np.isfinite(cl)) or np.any(cl < 0):
        raise ValueError('Invalid power spectrum')
    theta_grid = np.linspace(0, np.pi, args.theta_samples)
    xp, xm = xi_from_cl(theta_grid, cl)
    # Quantify interpolation error at midpoints against direct harmonic sums.
    rng = np.random.default_rng(20261002)
    check_theta = np.concatenate((rng.uniform(0, np.deg2rad(15), 300), rng.uniform(0, np.pi, 100)))
    cp, cm = xi_from_cl(check_theta, cl)
    interpolation_error = max(np.max(abs(cp - np.interp(check_theta,theta_grid,xp))),
                              np.max(abs(cm - np.interp(check_theta,theta_grid,xm)))) / xp[0]
    if interpolation_error > 1e-4:
        raise ValueError(f'Interpolation error {interpolation_error}; increase --theta-samples')
    print(f'pixels={len(pix)}, lmax={lmax}, xi_plus_zero={xp[0]:.10g}, interpolation/xi0={interpolation_error:.3g}',flush=True)
    signal = signal_covariance(vectors, theta_grid, xp, xm, xp[0])
    if not np.allclose(signal, signal.transpose(1,0,3,2), atol=xp[0]*1e-12, rtol=1e-12):
        raise AssertionError('Signal covariance is not transpose symmetric')
    counts, sums, paths, hashes, means, boundary_counts = [], [], [], [], [], []
    for b in range(len(edges)-1):
        path = args.pair_base/f'pairs_merged/pairs_bin_{b:02d}.npy'
        pairs = np.load(path, mmap_mode='r')
        if pairs.ndim != 2 or pairs.shape[1] != 2 or len(pairs)==0:
            raise ValueError(f'Invalid/empty pair file: {path}')
        i,j = pairs.T
        if np.any(i<0) or np.any(j>=len(pix)) or np.any(i>=j):
            raise ValueError(f'Pairs must be valid unordered i<j: {path}')
        keys = i.astype(np.int64)*len(pix)+j
        if len(np.unique(keys)) != len(keys):
            raise ValueError(f'Duplicate pairs: {path}')
        q=weights[i]*weights[j]
        angle=np.arccos(np.clip(np.sum(vectors[i]*vectors[j],axis=1),-1,1))
        degree=np.rad2deg(angle)
        boundary_counts.append(np.count_nonzero((degree<edges[b])|(degree>=edges[b+1])))
        geo=pair_geometry(pairs,vectors,weights)
        hij=np.einsum('nac,ncd,nbd->nab',geo[2],signal[i,j],geo[3],optimize=True)
        means.append(np.sum(q[:,None]*np.stack((hij[:,0,0]+hij[:,1,1],hij[:,0,0]-hij[:,1,1]),axis=1),axis=0)/q.sum())
        counts.append(len(pairs)); sums.append(q.sum()); paths.append(str(path.resolve()));hashes.append(file_hash(path))
        print(f'bin={b} full_pairs={len(pairs)} old-list boundary differences={boundary_counts[-1]}',flush=True)
    exclusive_save(args.output, vectors=vectors, pix_idx=pix, weights=weights,
        ell=np.arange(len(cl)),cl_bb=np.zeros_like(cl),
        nside=nside,zbin=zbin,theta=np.sqrt(edges[:-1]*edges[1:]),theta_edges=edges,
        theta_grid=theta_grid,xi_plus=xp,xi_minus=xm,xi_plus_zero=xp[0],cl_ee=cl,
        signal=signal,full_weight_sums=np.array(sums),full_counts=np.array(counts),
        pair_paths=np.array(paths),pair_hashes=np.array(hashes),
        metadata_sha256=file_hash(meta_path),cl_source=str(args.cl_file.resolve()),
        cl_source_sha256=file_hash(args.cl_file),cl_key=args.cl_key,
        prediction_mean=np.array(means).T.reshape(-1), boundary_counts=boundary_counts,
        interpolation_error_over_xi0=interpolation_error,
        convention='TreeCorr west/north, gamma=(healpy Q,-healpy U)',
        spectrum_note='Input numeric spectrum treated as EE directly, BB=0; pixwin=False; fwhm=0; no kappa-to-shear conversion',
        estimator='fixed original convergence full pixel-pair lists and weights; no shape noise')
    print(f'Saved {args.output}',flush=True)


if __name__=='__main__': main()
