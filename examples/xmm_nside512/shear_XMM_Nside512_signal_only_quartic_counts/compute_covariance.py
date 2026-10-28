#!/usr/bin/env python3
"""Full pair matrices or exact HT p=1/1000; four shear blocks computed together."""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from shear_core import (file_hash, exclusive_save, checkpoint_save, pair_geometry,
                        bernoulli_pairs, identical_sum, raw_pair_sum, assemble_blocks)
ROOT=Path(__file__).resolve().parent


def load_pairs(model, override=None):
    pairs=[]
    for b,(old,expected) in enumerate(zip(model['pair_paths'],model['pair_hashes'])):
        path=Path(old) if override is None else Path(override)/f'pairs_bin_{b:02d}.npy'
        if file_hash(path)!=expected:
            raise ValueError(f'Pair content differs from prepared model: {path}')
        pairs.append(np.load(path,mmap_mode='r'))
    return pairs


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--model',type=Path,default=ROOT/'outputs/model.npz')
    ap.add_argument('--pairs-dir',type=Path)
    ap.add_argument('--mode',choices=['full','thinned'],required=True)
    ap.add_argument('--row',type=int,help='Full mode: compute one angular-bin row (upper triangle)')
    ap.add_argument('--realizations',type=int,default=100)
    ap.add_argument('--base-seed',type=int,default=20260828)
    ap.add_argument('--left-chunk',type=int,default=0)
    ap.add_argument('--right-chunk',type=int,default=2000)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    if args.realizations<1 or args.left_chunk<0 or args.right_chunk<1:
        raise ValueError('Invalid realization/chunk parameter')
    if args.row is not None and args.mode!='full': raise ValueError('--row requires full')
    tag=('full' if args.row is None else f'full_row_{args.row:02d}') if args.mode=='full' else 'ht_p0.001'
    output_root = ROOT if args.mode=='full' else ROOT.parent/'shear_XMM_Nside512_signal_only_quartic_counts_down_sampling'
    output=args.output or output_root/f'outputs/{tag}.npz'
    if output.exists(): raise FileExistsError(f'Refusing existing final output: {output}')
    model_hash=file_hash(args.model)
    with np.load(args.model) as d: model={k:d[k] for k in d.files}
    pairs=load_pairs(model,args.pairs_dir)
    n=len(pairs)
    if args.row is not None and not 0<=args.row<n: raise ValueError('Invalid row')
    p=1. if args.mode=='full' else .001
    repetitions=1 if p==1 else args.realizations
    config=json.dumps(dict(model_hash=model_hash,mode=args.mode,row=args.row,n_realizations=repetitions,
                           seed=args.base_seed,p=p),sort_keys=True)
    cp=output.with_suffix('.checkpoint.npz')
    completed=0
    all_blocks=np.full((repetitions,n,n,2,2),np.nan)
    retained=np.zeros((repetitions,n),int)
    seeds=np.zeros((repetitions,n),np.uint64)
    if cp.exists():
        with np.load(cp) as d:
            if str(d['config'])!=config: raise ValueError('Checkpoint settings/model mismatch')
            completed=int(d['completed']); all_blocks[:]=d['blocks'];retained[:]=d['retained'];seeds[:]=d['seeds']
        print(f'Resuming {completed}/{repetitions}',flush=True)
    sequences=np.random.SeedSequence(args.base_seed).spawn(repetitions)
    started=time.monotonic()
    for r in range(completed,repetitions):
        selected=[]
        for b,seq in enumerate(sequences[r].spawn(n)):
            seeds[r,b]=seq.generate_state(1,dtype=np.uint64)[0]
            selected.append(bernoulli_pairs(pairs[b],p,np.random.default_rng(seq)))
            retained[r,b]=len(selected[-1])
        geometry=[pair_geometry(x,model['vectors'],model['weights']) for x in selected]
        diag=[identical_sum(g,model['signal']) for g in geometry] if p<1 else None
        rows=range(n) if args.row is None else [args.row]
        for a in rows:
            for b in range(a,n):
                if np.all(np.isfinite(all_blocks[r,a,b])): continue
                print(f'realization={r+1}/{repetitions} bin=({a},{b}) sizes={len(selected[a])}x{len(selected[b])}',flush=True)
                raw=raw_pair_sum(geometry[a],geometry[b],model['signal'],args.right_chunk,args.left_chunk)
                corrected=raw if p==1 else (raw-(1-p)*diag[a] if a==b else raw)/(p*p)
                block=corrected/(model['full_weight_sums'][a]*model['full_weight_sums'][b])
                all_blocks[r,a,b]=block
                all_blocks[r,b,a]=block.T  # +- is NOT symmetric on its own
                # Cell checkpoints also allow long full rows to resume.
                checkpoint_save(cp,config=config,completed=r,blocks=all_blocks,retained=retained,seeds=seeds)
                print(f'cell_done elapsed={time.monotonic()-started:.1f}s',flush=True)
        checkpoint_save(cp,config=config,completed=r+1,blocks=all_blocks,retained=retained,seeds=seeds)
    mean=np.mean(all_blocks,axis=0)
    sem=np.std(all_blocks,axis=0,ddof=1)/np.sqrt(repetitions) if repetitions>1 else np.zeros_like(mean)
    joint=assemble_blocks(mean)
    if args.row is None and (not np.all(np.isfinite(joint)) or not np.allclose(joint,joint.T,atol=1e-20,rtol=1e-10)):
        raise AssertionError('Final covariance invalid')
    exclusive_save(output, covariance_joint=joint,covariance_sem=assemble_blocks(sem),
        covariance_realizations=np.array([assemble_blocks(x) for x in all_blocks]),
        blocks=mean,retained_counts=retained,bin_seeds=seeds,keep_probability=p,
        n_realizations=repetitions,base_seed=args.base_seed,row=-1 if args.row is None else args.row,
        model_sha256=model_hash,config=config,theta=model['theta'],theta_edges=model['theta_edges'],
        full_weight_sums=model['full_weight_sums'],full_counts=model['full_counts'],
        ordering='plus bins, then minus bins',estimator=model['estimator'])
    print(f'Saved {output}',flush=True)


if __name__=='__main__':main()
