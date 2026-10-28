#!/usr/bin/env python3
"""Read-only server checks before preparing/submitting; never runs qsub."""
import argparse
import importlib
import os
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parent
SIM_ROOT=ROOT.parent/'shear_XMM_Nside512_signal_only_sims'
THIN_ROOT=ROOT.parent/'shear_XMM_Nside512_signal_only_quartic_counts_down_sampling'
from compute_covariance import load_pairs
from shear_core import file_hash

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--stage',choices=['prepare','thinned','simulation'],required=True)
    ap.add_argument('--pair-base',type=Path,default=ROOT/'outputs/precompute')
    ap.add_argument('--cl-file',type=Path,
                    help='Required for prepare; input spectrum NPZ')
    ap.add_argument('--model',type=Path,default=ROOT/'outputs/model.npz')
    ap.add_argument('--catalog',type=Path,
                    help='Required for simulation; authorized catalog FITS')
    ap.add_argument('--grid-only',action='store_true')
    ap.add_argument('--require-conda',action='store_true')
    args=ap.parse_args()
    if args.require_conda and os.environ.get('CONDA_DEFAULT_ENV')!='cellcov_3':
        raise RuntimeError('Activate conda cellcov_3 before using PBS scripts')
    for name in ['numpy','healpy','treecorr','astropy']:
        mod=importlib.import_module(name)
        print(f'{name}={getattr(mod,"__version__","unknown")}')
    print(f'python={sys.executable}; conda={os.environ.get("CONDA_DEFAULT_ENV","unset")}')
    print('Expected PBS: prepare 1 CPU/16GB/3h; HT thinned 1 CPU/16GB/24h; simulation 8 CPU/32GB/72h')
    print('Queue/resource availability and installed server package versions must be checked on the server.')
    if args.stage=='prepare':
        if args.cl_file is None:
            raise ValueError('--cl-file is required for prepare')
        paths=[args.pair_base/'pair_meta.npz',args.cl_file]+[args.pair_base/f'pairs_merged/pairs_bin_{b:02d}.npy' for b in range(10)]
        for p in paths:
            if not p.is_file(): raise FileNotFoundError(p)
        if args.model.exists(): raise FileExistsError(f'Model already prepared: {args.model}')
    else:
        with np.load(args.model) as d:m={k:d[k] for k in d.files}
        pairs=load_pairs(m)
        print('model_sha256=',file_hash(args.model));print('full_counts=',m['full_counts'])
        if int(m['nside'])!=512 or int(m['zbin'])!=4:raise ValueError('Production model must be Nside512/bin4')
        if args.stage=='simulation':
            if not args.grid_only and args.catalog is None:
                raise ValueError('--catalog is required for simulation')
            if not args.grid_only and not args.catalog.is_file():raise FileNotFoundError(args.catalog)
        if args.stage=='thinned' and (THIN_ROOT/'outputs/ht_p0.001.npz').exists():raise FileExistsError('HT final output already exists')
        if args.stage=='simulation' and (SIM_ROOT/'outputs/shear_signal_3000sims.npz').exists():raise FileExistsError('Simulation final output already exists')
    print('Preflight passed (no job submitted).')
if __name__=='__main__':main()
