#!/usr/bin/env python3
"""Small synthetic catalog integration test; no production files are written."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import hashlib
import numpy as np
import healpy as hp
from astropy.io import fits

ROOT = Path(__file__).resolve().parent
LENSING = ROOT.parent
spec = importlib.util.spec_from_file_location('noise_sim', ROOT/'simulate_shear_shape_noise.py')
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)


def run(args, success=True):
    result = subprocess.run([sys.executable, *map(str,args)], capture_output=True, text=True)
    if (result.returncode == 0) != success:
        raise AssertionError(result.stdout + result.stderr)
    return result


def main():
    with tempfile.TemporaryDirectory(prefix='shear_noise_test_') as temp:
        d = Path(temp); rng = np.random.default_rng(51); ngal=400
        ra=34+rng.uniform(0,4,ngal); dec=-5+rng.uniform(0,4,ngal)
        w=rng.uniform(.5,2,ngal); e1=rng.normal(0,.2,ngal); e2=rng.normal(0,.25,ngal)
        catalog=d/'catalog.fits'
        cols=[fits.Column(name=k,format=fmt,array=v) for k,fmt,v in [
            ('i_ra','D',ra),('i_dec','D',dec),('i_hsmshaperegauss_derived_weight','D',w),
            ('e1','D',e1),('e2','D',e2),('b_mode_mask','L',np.ones(ngal,bool)),
            ('hsc_y3_zbin','K',np.full(ngal,4))]]
        fits.BinTableHDU.from_columns(cols).writeto(catalog)
        pix=hp.ang2pix(512,ra,dec,lonlat=True); ids=np.unique(pix)
        weights=np.bincount(pix,weights=w,minlength=hp.nside2npix(512))[ids]
        vec=np.array(hp.pix2vec(512,ids)).T
        i,j=np.triu_indices(len(ids),1); angle=np.rad2deg(np.arccos(np.clip(np.sum(vec[i]*vec[j],axis=1),-1,1)))
        edges=np.geomspace(.1,5,11); paths=[]; hashes=[]; sums=[]
        for b in range(10):
            keep=(angle>=edges[b])&(angle<edges[b+1]); pairs=np.c_[i[keep],j[keep]]
            assert len(pairs)>0
            path=d/f'pairs_{b}.npy';np.save(path,pairs)
            paths.append(str(path));hashes.append(hashlib.sha256(path.read_bytes()).hexdigest())
            sums.append(np.sum(weights[pairs[:,0]]*weights[pairs[:,1]]))
        cl=np.zeros(1536);cl[2:]=1e-8/(np.arange(2,1536)+1)**1.5
        model=d/'model.npz'
        np.savez(model,cl_ee=cl,cl_bb=np.zeros_like(cl),ell=np.arange(1536),nside=512,zbin=4,
                 theta=np.sqrt(edges[:-1]*edges[1:]),theta_edges=edges,pix_idx=ids,weights=weights,
                 vectors=vec,pair_paths=paths,pair_hashes=hashes,full_weight_sums=sums,estimator='test')
        # Run the existing, unchanged signal-only program as the regression reference.
        signal=d/'signal.npz'
        run([LENSING/'shear_XMM_Nside512_signal_only_sims/simulate_shear.py',
             '--model',model,'--catalog',catalog,'--nsims','2','--threads','1','--output',signal])
        # Independently reproduce the Convergence second-moment formula.
        r1=1-np.sum(w*e1**2)/w.sum();r2=1-np.sum(w*e2**2)/w.sum()
        s1=np.sum(w*(e1/(2*r1))**2)/w.sum();s2=np.sum(w*(e2/(2*r2))**2)/w.sum()
        conv=d/'convergence.npz'
        np.savez(conv,n_galaxies=ngal,nside=512,zbin=4,responsivity_1=r1,responsivity_2=r2,
                 sigma_g1_squared=s1,sigma_g2_squared=s2,sigma_kappa_noise_squared=(s1+s2)/2)
        out=d/'total.npz'
        args=[ROOT/'simulate_shear_shape_noise.py','--model',model,'--catalog',catalog,
              '--reference-signal',signal,'--reference-convergence',conv,'--nsims','2','--threads','1']
        run([*args,'--preflight'])
        run([*args,'--output',out])
        with np.load(out) as f, np.load(signal) as old:
            np.testing.assert_allclose(f['xi_signal_sims'],old['xi_catalog'],rtol=1e-9,atol=1e-14)
            for label in mod.LABELS:
                assert f[f'xi_{label}_sims'].shape==(2,20)
                np.testing.assert_allclose(f[f'covariance_{label}'],np.cov(f[f'xi_{label}_sims'],rowvar=False))
            baseline={k:f[k] for k in f.files}
        # Historical outputs are checks, not inputs to the stochastic calculation.
        independent=d/'independent.npz'
        result=run([ROOT/'simulate_shear_shape_noise.py','--model',model,'--catalog',catalog,
                    '--nsims','2','--threads','1','--output',independent])
        assert 'signal_reference=NOT_REQUESTED' in result.stdout
        assert 'Convergence reference: NOT REQUESTED' in result.stdout
        with np.load(independent) as f:
            import json
            cfg=json.loads(str(f['config']))
            assert cfg['signal_reference_sha256'] is None
            assert cfg['convergence_reference_sha256'] is None
            for key, value in baseline.items():
                if key != 'config':
                    np.testing.assert_array_equal(f[key], value, err_msg=key)
        resumed=d/'resumed.npz'
        with np.load(out.with_suffix('.checkpoint.npz')) as f: state={k:f[k] for k in f.files}
        state['completed']=1
        for k in state:
            if k.startswith('xi_'): state[k][1:]=np.nan
        np.savez(resumed.with_suffix('.checkpoint.npz'),**state)
        run([*args,'--output',resumed])
        with np.load(resumed) as f:
            for label in mod.LABELS:
                np.testing.assert_array_equal(f[f'xi_{label}_sims'],baseline[f'xi_{label}_sims'])
        run([*args,'--output',out],success=False)
        run([*args,'--seed-start','100','--preflight'],success=False)
        np.random.seed(9); before=np.random.get_state(); noise=mod.draw_noise(300000,(s1+s2)/2,0)
        after=np.random.get_state();np.testing.assert_array_equal(before[1],after[1]);assert before[2:]==after[2:]
        np.testing.assert_allclose(noise.var(axis=1),(s1+s2)/2,rtol=.015)
        assert abs(np.corrcoef(noise)[0,1])<.01
        assert not np.array_equal(noise[:,0],noise[:,1])
        np.testing.assert_array_equal(noise,mod.draw_noise(300000,(s1+s2)/2,0))
        print('PASS: existing signal regression; noise convention; 3 catalog outputs; covariance; identical outputs without historical references; resume; overwrite/seed guards; independent RNG.')


def test_catalog_simulation():
    main()


if __name__=='__main__':
    main()
