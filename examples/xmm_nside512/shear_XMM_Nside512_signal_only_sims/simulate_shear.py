#!/usr/bin/env python3
"""Nside512 signal-only maps: exact saved-grid estimator and catalog TreeCorr."""
import argparse
import json
import time
from pathlib import Path
import numpy as np
import healpy as hp
import treecorr
from astropy.io import fits
import sys
COMMON=Path(__file__).resolve().parents[1]/'shear_XMM_Nside512_signal_only_quartic_counts'
sys.path.insert(0,str(COMMON))
from shear_core import file_hash, exclusive_save,checkpoint_save,pair_geometry,measure_grid
from compute_covariance import load_pairs
ROOT=Path(__file__).resolve().parent
CATALOG=Path(__file__).resolve().parents[1]/'data'/'XMM.fits'


def load_catalog(path,nside,zbin):
    with fits.open(path,memmap=True) as hdul:
        t=hdul[1].data
        sel=(t['b_mode_mask']==True)&(t['hsc_y3_zbin']==zbin)
        ra=np.asarray(t['i_ra'][sel],float);dec=np.asarray(t['i_dec'][sel],float)
        w=np.asarray(t['i_hsmshaperegauss_derived_weight'][sel],float)
    good=np.isfinite(ra)&np.isfinite(dec)&np.isfinite(w)&(w>0)
    ra,dec,w=ra[good],dec[good],w[good]
    pix=hp.ang2pix(nside,ra,dec,lonlat=True)
    return ra,dec,w,pix


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--model',type=Path,default=COMMON/'outputs/model.npz')
    ap.add_argument('--pairs-dir',type=Path)
    ap.add_argument('--catalog',type=Path,default=CATALOG)
    ap.add_argument('--grid-only',action='store_true')
    ap.add_argument('--nsims',type=int,default=3000)
    ap.add_argument('--seed-start',type=int,default=0)
    ap.add_argument('--checkpoint-every',type=int,default=1)
    ap.add_argument('--threads',type=int,default=1)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    if args.nsims<2 or args.seed_start<0 or args.threads<1 or args.checkpoint_every<1:
        raise ValueError('Invalid simulation parameters')
    output=args.output or ROOT/f'outputs/shear_signal_{args.nsims}sims.npz'
    if output.exists(): raise FileExistsError(output)
    with np.load(args.model) as f: m={k:f[k] for k in f.files}
    pairs=load_pairs(m,args.pairs_dir)
    geometry=[pair_geometry(p,m['vectors'],m['weights']) for p in pairs]
    n=len(pairs); nside=int(m['nside']);cl=m['cl_ee'];zero=np.zeros_like(cl)
    if cl.ndim != 1 or len(cl) < 3 or np.any(cl[:2] != 0):
        raise ValueError('Model EE must be indexed by ell with C_0=C_1=0')
    if 'ell' in m and not np.array_equal(m['ell'],np.arange(len(cl))):
        raise ValueError('Model multipole indices are inconsistent')
    if 'cl_bb' in m and (m['cl_bb'].shape != cl.shape or np.any(m['cl_bb'] != 0)):
        raise ValueError('Signal-only E-mode model requires BB=0')
    catalog_hash='none'; ngal=0
    if not args.grid_only:
        ra,dec,w,pix=load_catalog(args.catalog,nside,int(m['zbin']))
        ngal=len(w);catalog_hash=file_hash(args.catalog)
        weight_map=np.bincount(pix,weights=w,minlength=hp.nside2npix(nside))
        if np.count_nonzero(weight_map)!=len(m['pix_idx']) or not np.allclose(weight_map[m['pix_idx']],m['weights'],rtol=3e-6,atol=1e-5):
            raise ValueError('Catalog selection does not reproduce grid weights')
        # Preserve the historical nearest-pixel component assignment explicitly.
        # No transport from pixel centre to galaxy: this is a separate diagnostic estimator.
        edges=m['theta_edges']
        template=treecorr.GGCorrelation(min_sep=edges[0],max_sep=edges[-1],nbins=n,
            sep_units='deg')
    config=json.dumps(dict(model_sha256=file_hash(args.model),catalog_sha256=catalog_hash,
        nsims=args.nsims,seed_start=args.seed_start,grid_only=args.grid_only,
        measurement_policy='treecorr_defaults_v2',treecorr_version=treecorr.__version__,
        metric='Euclidean',bin_slop=None if args.grid_only else template.bin_slop,
        angle_slop=None if args.grid_only else template.angle_slop,pixwin=False,fwhm=0),sort_keys=True)
    grid=np.full((args.nsims,2*n),np.nan);cat=np.full_like(grid,np.nan)
    imag=np.full_like(grid,np.nan);completed=0
    cp=output.with_suffix('.checkpoint.npz')
    if cp.exists():
        with np.load(cp) as d:
            if str(d['config'])!=config: raise ValueError('Checkpoint configuration mismatch')
            grid[:]=d['grid'];cat[:]=d['catalog'];imag[:]=d['imag'];completed=int(d['completed'])
    print(f'measurement config={config}; resumed={completed}',flush=True)
    started=time.monotonic()
    for r in range(completed,args.nsims):
        tick=time.monotonic()
        seed=args.seed_start+r
        np.random.seed(seed)
        maps=hp.synfast([zero,cl,zero,zero],nside,lmax=len(cl)-1,new=True,pol=True,pixwin=False,fwhm=0.)
        map_done=time.monotonic()
        g=np.column_stack((maps[1][m['pix_idx']],-maps[2][m['pix_idx']]))
        grid[r]=measure_grid(g,geometry,m['full_weight_sums'])
        grid_done=time.monotonic()
        if not args.grid_only:
            gg=template.copy()
            gg.process(treecorr.Catalog(ra=ra,dec=dec,ra_units='deg',dec_units='deg',
                w=w,g1=maps[1][pix],g2=-maps[2][pix]),num_threads=args.threads)
            cat[r]=np.concatenate((gg.xip,gg.xim));imag[r]=np.concatenate((gg.xip_im,gg.xim_im))
        catalog_done=time.monotonic()
        print(f'timing map={map_done-tick:.2f}s grid={grid_done-map_done:.2f}s catalog={catalog_done-grid_done:.2f}s',flush=True)
        done=r+1
        if done%args.checkpoint_every==0 or done==args.nsims:
            checkpoint_save(cp,config=config,completed=done,grid=grid,catalog=cat,imag=imag)
        elapsed=time.monotonic()-started
        print(f'realization={done}/{args.nsims} seed={seed} elapsed={elapsed:.1f}s ETA={elapsed/(done-completed)*(args.nsims-done):.1f}s',flush=True)
    extra={}
    if not args.grid_only:
        extra=dict(xi_catalog=cat,xi_catalog_imag=imag,covariance_catalog=np.cov(cat,rowvar=False,ddof=1),mean_catalog=cat.mean(axis=0))
    exclusive_save(output,xi_grid=grid,covariance_grid=np.cov(grid,rowvar=False,ddof=1),mean_grid=grid.mean(axis=0),
        seeds=np.arange(args.seed_start,args.seed_start+args.nsims),nsims=args.nsims,
        model_sha256=file_hash(args.model),config=config,theta=m['theta'],theta_edges=m['theta_edges'],
        nside=nside,zbin=m['zbin'],cl_ee=cl,n_galaxies=ngal,ordering='plus bins, then minus bins',
        grid_estimator=m['estimator'],catalog_estimator='galaxy-position TreeCorr defaults (Euclidean, default slops); nearest-pixel Q,-U without component transport; distinct from grid',
        **extra)
    print(f'Saved {output}',flush=True)


if __name__=='__main__':main()
