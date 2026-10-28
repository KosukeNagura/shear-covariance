#!/usr/bin/env python3
"""Rebuild pairs/model/simulations with packaged scientific scripts; compare separately.

Historical arrays are read only for comparison. Only bin definitions and the
interpolation-grid size are extracted as settings from the historical model.
No historical pairs, model covariance, or realizations are generation inputs.
"""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

import numpy as np

PACKAGE = Path(__file__).resolve().parents[1]
EXAMPLES = PACKAGE / 'examples/xmm_nside512'
COMMON = EXAMPLES / 'shear_XMM_Nside512_signal_only_quartic_counts'
PAIR = EXAMPLES / 'pair_precompute_existing'
SIGNAL = EXAMPLES / 'shear_XMM_Nside512_signal_only_sims/simulate_shear.py'
NOISE = EXAMPLES / 'shear_XMM_Nside512_signal_plus_shape_noise_sims/simulate_shear_shape_noise.py'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def compare_array(actual, expected, *, rtol=0., atol=0.):
    """Always report exact equality, also when accepting floating tolerance."""
    a, b = np.asarray(actual), np.asarray(expected)
    result = dict(shape=list(a.shape), reference_shape=list(b.shape),
                  rtol=rtol, atol=atol, exact=False, passed=False)
    if a.shape != b.shape:
        return result
    if a.dtype.kind in 'biufc' and b.dtype.kind in 'biufc':
        if not (np.isfinite(a).all() and np.isfinite(b).all()):
            result['nonfinite'] = True
            return result
        result['exact'] = bool(np.array_equal(a, b))
        if a.dtype.kind in 'biu' and b.dtype.kind in 'biu':
            result['passed'] = result['exact']
            return result
        result['max_abs_difference'] = float(np.max(np.abs(a.astype(float)-b.astype(float)))) if a.size else 0.
        result['passed'] = bool(np.allclose(a, b, rtol=rtol, atol=atol))
    else:
        result['passed'] = result['exact'] = bool(np.array_equal(a, b))
    return result


def matched_rows(new_seeds, old_seeds):
    old = np.asarray(old_seeds)
    if len(np.unique(old)) != len(old):
        raise ValueError('Duplicate reference seeds')
    rows = {int(seed): i for i, seed in enumerate(old)}
    return [rows[int(seed)] for seed in new_seeds]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for option in ('catalog', 'cl-file', 'reference-model', 'reference-signal', 'reference-noise', 'work-root'):
        ap.add_argument('--'+option, required=True, type=Path)
    ap.add_argument('--cl-key', default='cl_theory')
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--nsims', type=int, default=3000)
    args = ap.parse_args()
    if args.nsims < 3 or args.threads < 1:
        raise ValueError('Require at least 3 simulations and one thread')
    work = args.work_root.resolve()
    # Never mix a previous partial run with a fresh end-to-end validation.
    work.mkdir(parents=True, exist_ok=False)
    report = dict(status='running', stages=[], comparisons=[],
                  scope='catalog -> regenerated pairs -> regenerated model -> regenerated simulations',
                  pending=['thinned covariance from regenerated model/simulations', 'notebook using regenerated outputs'])
    report_path = work/'validation_report.json'

    def save():
        temp = report_path.with_suffix('.tmp')
        temp.write_text(json.dumps(report, indent=2, ensure_ascii=False)+'\n')
        temp.replace(report_path)

    def check(label, actual, expected, **tolerance):
        result = dict(name=label, **compare_array(actual, expected, **tolerance))
        report['comparisons'].append(result)
        save()
        print(('PASS' if result['passed'] else 'FAIL') + ': ' + label +
              ' exact=' + str(result['exact']), flush=True)
        if not result['passed']:
            raise AssertionError(f'{label}: {result}')

    env = dict(os.environ)
    # Do not accidentally import an external 2PCF or a user PYTHONPATH copy.
    env.pop('PYTHONPATH', None)
    env['PYTHONNOUSERSITE'] = '1'
    env['SHEAR_2PCF_FUNCS'] = str(PAIR/'legacy_2pcf_funcs')
    env['SHEAR_PAIR_WORK_ROOT'] = str(work/'pairs')
    env['SHEAR_PAIR_CATALOG'] = str(args.catalog.resolve())
    env['SHEAR_PAIR_THETA_NPZ'] = str(work/'bin_settings.npz')

    def run(label, script, *options):
        command = [sys.executable, '-u', str(script), *map(str, options)]
        stage = dict(name=label, command=command, status='running')
        report['stages'].append(stage)
        save()
        print('\nSTART '+label+'\n'+shlex.join(command), flush=True)
        start = time.monotonic()
        proc = subprocess.run(command, cwd=work, env=env, check=False)
        stage.update(exit_status=proc.returncode, seconds=time.monotonic()-start,
                     status='passed' if proc.returncode == 0 else 'failed')
        save()
        if proc.returncode:
            raise RuntimeError(f'{label} exited with {proc.returncode}')

    def compare_simulations(signal_path, noise_path, label):
        with np.load(signal_path) as new, np.load(args.reference_signal) as old:
            rows = matched_rows(new['seeds'], old['seeds'])
            nc, oc = json.loads(str(new['config'])), json.loads(str(old['config']))
            for key in ('catalog_sha256', 'metric', 'bin_slop', 'angle_slop', 'treecorr_version', 'pixwin', 'fwhm'):
                check(label+'/signal/config/'+key, nc[key], oc[key])
            for key in ('theta', 'theta_edges', 'cl_ee', 'nside', 'zbin', 'n_galaxies'):
                check(label+'/signal/'+key, new[key], old[key])
            for key in ('xi_grid', 'xi_catalog', 'xi_catalog_imag'):
                check(label+'/signal/'+key, new[key], old[key][rows], rtol=1e-9, atol=1e-14)
            check(label+'/signal/covariance_catalog', new['covariance_catalog'],
                  np.cov(old['xi_catalog'][rows], rowvar=False, ddof=1), rtol=2e-9, atol=1e-20)
        with np.load(noise_path) as new, np.load(args.reference_noise) as old:
            rows = matched_rows(new['seeds'], old['seeds'])
            nc, oc = json.loads(str(new['config'])), json.loads(str(old['config']))
            for key in ('catalog_sha256', 'metric', 'bin_slop', 'angle_slop', 'treecorr_version',
                        'numpy_version', 'healpy_version', 'pixwin', 'fwhm', 'noise_seed_offset'):
                check(label+'/noise/config/'+key, nc[key], oc[key])
            for key in ('theta', 'theta_edges', 'cl_ee', 'nside', 'zbin', 'n_galaxies', 'noise_seed_offset'):
                check(label+'/noise/'+key, new[key], old[key])
            for key in ('responsivity_1', 'responsivity_2', 'sigma_g1_squared',
                        'sigma_g2_squared', 'sigma_component_squared'):
                check(label+'/noise/'+key, new[key], old[key], rtol=1e-12)
            for component in ('signal', 'noise', 'total'):
                for suffix in ('', '_imag'):
                    key = f'xi_{component}{suffix}_sims'
                    check(label+'/noise/'+key, new[key], old[key][rows], rtol=1e-9, atol=1e-14)
                check(label+'/noise/covariance_'+component, new['covariance_'+component],
                      np.cov(old[f'xi_{component}_sims'][rows], rowvar=False, ddof=1),
                      rtol=2e-9, atol=1e-20)

    try:
        inputs = {name: getattr(args, name).resolve() for name in
                  ('catalog', 'cl_file', 'reference_model', 'reference_signal', 'reference_noise')}
        report['inputs'] = {name: dict(path=str(path), sha256=digest(path)) for name, path in inputs.items()}
        report['python'] = sys.version
        report['packages'] = {name: importlib.metadata.version(name) for name in
                              ('numpy', 'scipy', 'healpy', 'astropy', 'treecorr')}
        report['code_sha256'] = {str(p.relative_to(PACKAGE)): digest(p)
                                 for p in sorted(PACKAGE.rglob('*.py'))}
        report['settings'] = dict(nsims=args.nsims, threads=args.threads, seed_start=0,
                                  pair_chunk_size=2000, cl_key=args.cl_key)
        with np.load(args.reference_model) as old:
            # ONLY analysis settings are exported. No weights, pairs, xi or S.
            np.savez(work/'bin_settings.npz', theta=old['theta'], theta_edges=old['theta_edges'])
            samples = len(old['theta_grid'])
            check('input spectrum checksum', report['inputs']['cl_file']['sha256'], str(old['cl_source_sha256']))
            check('input spectrum key', args.cl_key, str(old['cl_key']))
        report['settings']['theta_samples'] = samples
        with np.load(args.reference_noise) as old:
            check('catalog checksum', report['inputs']['catalog']['sha256'], str(old['catalog_sha256']))
        save()
        run('synthetic simulation regression including no-reference mode', NOISE.with_name('test_catalog_simulation.py'))
        run('catalog to pixel metadata', PAIR/'prep_pair_metadata.py')
        precompute = work/'pairs/outputs/precompute'
        with np.load(precompute/'pair_meta.npz') as new, np.load(args.reference_model) as old:
            for nk, ok in [('pix_idx','pix_idx'), ('weights','weights'), ('theta','theta'),
                           ('theta_edges','theta_edges'), ('nside','nside'), ('iz','zbin')]:
                check('metadata/'+nk, new[nk], old[ok])
            npixels = len(new['pix_idx'])
        for chunk_id, start in enumerate(range(0, npixels, 2000)):
            run(f'pair chunk {chunk_id}', PAIR/'pair_chunk_worker.py', chunk_id, start, min(start+2000, npixels))
            for b in range(10):
                if not (precompute/f'pair_chunks/pairs_bin_{b:02d}_chunk_{chunk_id:04d}.npy').is_file():
                    raise FileNotFoundError(f'Missing chunk {chunk_id}, bin {b}')
        run('merge complete pair chunks', PAIR/'merge_pair_chunks.py')
        with np.load(args.reference_model) as old:
            for b, (path, expected_hash) in enumerate(zip(old['pair_paths'], old['pair_hashes'])):
                check(f'reference pair integrity/{b}', digest(str(path)), str(expected_hash))
                check(f'pair arrays/{b}', np.load(precompute/f'pairs_merged/pairs_bin_{b:02d}.npy'),
                      np.load(str(path)))
            check('pair counts', np.load(precompute/'pair_counts.npy'), old['full_counts'])
        model = work/'model.npz'
        run('model from regenerated pairs and input spectrum', COMMON/'prepare_shear.py',
            '--pair-base', precompute, '--cl-file', args.cl_file.resolve(), '--cl-key', args.cl_key,
            '--theta-samples', samples, '--output', model)
        with np.load(model) as new, np.load(args.reference_model) as old:
            # Absolute paths/archive hashes intentionally differ, numeric contents must reproduce.
            for key in ('pix_idx','nside','zbin','ell','theta','theta_edges','theta_grid','cl_ee','cl_bb',
                        'weights','vectors','full_counts','boundary_counts'):
                check('model/'+key, new[key], old[key])
            for key in ('xi_plus','xi_minus','xi_plus_zero','signal','full_weight_sums','prediction_mean',
                        'interpolation_error_over_xi0'):
                check('model/'+key, new[key], old[key], rtol=1e-12, atol=1e-18)
        for count, label in [(3, 'pilot')] + ([(args.nsims, 'production')] if args.nsims > 3 else []):
            signal_path = work/f'signal_{count}sims.npz'
            noise_path = work/f'noise_{count}sims.npz'
            run(label+' signal simulation', SIGNAL, '--model', model, '--catalog', args.catalog.resolve(),
                '--nsims', count, '--threads', args.threads, '--output', signal_path)
            # The noise run uses NO historical simulation files, including for preflight.
            run(label+' noise simulation', NOISE, '--model', model, '--catalog', args.catalog.resolve(),
                '--nsims', count, '--threads', args.threads, '--output', noise_path)
            compare_simulations(signal_path, noise_path, label)
        report['outputs'] = {str(p.relative_to(work)): digest(p) for p in work.glob('*.npz')}
        report['status'] = 'upstream_passed'
        save()
        print('PASS: regenerated metadata, pair arrays, model, signal/noise simulations and covariances.\n'
              'Still pending: thinned analytic covariance from these new inputs, then notebook.', flush=True)
    except Exception as error:
        report['status'] = 'failed'
        report['error'] = str(error)
        save()
        raise


if __name__ == '__main__':
    main()
