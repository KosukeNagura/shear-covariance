#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
from shear_core import exclusive_save,assemble_blocks
ROOT=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--directory',type=Path,default=ROOT/'outputs');a=ap.parse_args()
    output=a.directory/'full.npz'
    if output.exists():raise FileExistsError(output)
    with np.load(a.directory/'full_row_00.npz') as d:
        blocks=np.full_like(d['blocks'],np.nan); model=str(d['model_sha256'])
        theta=d['theta'];edges=d['theta_edges'];sums=d['full_weight_sums'];counts=d['full_counts']
    n=len(theta)
    for row in range(n):
        with np.load(a.directory/f'full_row_{row:02d}.npz') as d:
            if int(d['row'])!=row or str(d['model_sha256'])!=model or float(d['keep_probability'])!=1:
                raise ValueError('Row provenance mismatch')
            for b in range(row,n):
                blocks[row,b]=d['blocks'][row,b];blocks[b,row]=blocks[row,b].T
    joint=assemble_blocks(blocks)
    if not np.all(np.isfinite(joint)) or not np.allclose(joint,joint.T): raise ValueError('Invalid merge')
    exclusive_save(output,covariance_joint=joint,covariance_sem=np.zeros_like(joint),blocks=blocks,
        model_sha256=model,theta=theta,theta_edges=edges,full_weight_sums=sums,full_counts=counts,
        keep_probability=1.,n_realizations=1,ordering='plus bins, then minus bins')
    print(output)
if __name__=='__main__':main()
