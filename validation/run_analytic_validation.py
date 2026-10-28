#!/usr/bin/env python3
"""Validate both analytic HT covariances without waiting for any simulations."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

import numpy as np
from run_from_catalog import compare_array, digest

PACKAGE = Path(__file__).resolve().parents[1]
EXAMPLES = PACKAGE/'examples/xmm_nside512'
COMMON = EXAMPLES/'shear_XMM_Nside512_signal_only_quartic_counts'
NOISE = EXAMPLES/'shear_XMM_Nside512_signal_plus_shape_noise_quartic_counts/compute_covariance_noise.py'


def upstream_ready(path):
    """Require completed model/pair checks, not the ongoing simulation stages."""
    report = json.loads((path/'validation_report.json').read_text())
    checks = {c['name']: c for c in report['comparisons']}
    required = ['input spectrum checksum', 'input spectrum key', 'catalog checksum', 'pair counts']
    required += ['metadata/'+k for k in ('pix_idx', 'weights', 'theta', 'theta_edges', 'nside', 'iz')]
    required += [f'{prefix}/{b}' for prefix in ('reference pair integrity', 'pair arrays') for b in range(10)]
    required += ['model/'+k for k in ('pix_idx', 'nside', 'zbin', 'ell', 'theta', 'theta_edges',
                    'theta_grid', 'cl_ee', 'cl_bb', 'weights', 'vectors', 'full_counts', 'boundary_counts',
                    'xi_plus', 'xi_minus', 'xi_plus_zero', 'signal', 'full_weight_sums', 'prediction_mean',
                    'interpolation_error_over_xi0')]
    missing = [name for name in required if not checks.get(name, {}).get('passed', False)]
    if missing:
        raise ValueError('Upstream checks not passed: '+', '.join(missing))
    if not (path/'model.npz').is_file():
        raise FileNotFoundError(path/'model.npz')
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--upstream-run', required=True, type=Path)
    ap.add_argument('--reference-signal-ht', type=Path)
    ap.add_argument('--reference-noise-ht', type=Path)
    ap.add_argument('--work-root', type=Path)
    ap.add_argument('--check-ready', action='store_true', help='Only read upstream completion checks; no calculation')
    args = ap.parse_args()
    upstream = args.upstream_run.resolve()
    previous = upstream_ready(upstream)
    if args.check_ready:
        print('READY: catalog metadata, pair arrays and model checks passed; no simulation completion required.')
        return
    if any(x is None for x in (args.work_root, args.reference_signal_ht, args.reference_noise_ht)):
        ap.error('Supply --work-root and both reference HT paths')
    work = args.work_root.resolve()
    work.mkdir(parents=True, exist_ok=False)
    state = dict(status='running', scope='analytic covariance only; no simulation input',
                 comparisons=[], commands=[], pending=['full simulation validation and comparison plots'])

    def save():
        p = work/'analytic_report.json'
        tmp = p.with_suffix('.tmp')
        tmp.write_text(json.dumps(state, indent=2)+'\n')
        tmp.replace(p)

    def check(name, new, old, **tol):
        result = dict(name=name, **compare_array(new, old, **tol))
        state['comparisons'].append(result)
        save()
        print(f"{'PASS' if result['passed'] else 'FAIL'}: {name}; exact={result['exact']}", flush=True)
        if not result['passed']:
            raise AssertionError(result)

    def run(label, script, *options):
        cmd = [sys.executable, '-u', str(script), *map(str, options)]
        item = dict(name=label, command=cmd, status='running')
        state['commands'].append(item)
        save()
        print('\nSTART '+label+'\n'+shlex.join(cmd), flush=True)
        env = dict(os.environ, PYTHONNOUSERSITE='1')
        env.pop('PYTHONPATH', None)
        code = subprocess.run(cmd, cwd=work, env=env).returncode
        item.update(exit_status=code, status='passed' if code == 0 else 'failed')
        save()
        if code:
            raise RuntimeError(f'{label}: exit {code}')

    try:
        catalog = Path(previous['inputs']['catalog']['path'])
        model = upstream/'model.npz'
        reference_model = Path(previous['inputs']['reference_model']['path'])
        model_hash = digest(model)
        state['inputs'] = {name: dict(path=str(path), sha256=digest(path)) for name, path in (
            ('catalog', catalog), ('model', model), ('reference_model', reference_model),
            ('reference_signal_ht', args.reference_signal_ht), ('reference_noise_ht', args.reference_noise_ht))}
        state['code_sha256'] = {str(p.relative_to(PACKAGE)): digest(p) for p in sorted(PACKAGE.rglob('*.py'))}
        # Save the upstream evidence at this instant; never edit its live report.
        (work/'upstream_evidence.json').write_text(json.dumps(previous, indent=2)+'\n')
        check('catalog integrity', digest(catalog), previous['inputs']['catalog']['sha256'])
        check('historical model integrity', digest(reference_model), previous['inputs']['reference_model']['sha256'])
        with np.load(model) as fresh, np.load(reference_model) as old:
            for key in ('pix_idx','weights','vectors','cl_ee','cl_bb','theta','theta_edges','full_counts','full_weight_sums','signal'):
                check('model/'+key, fresh[key], old[key], rtol=1e-12, atol=1e-18)
            for b, (path, expected) in enumerate(zip(fresh['pair_paths'], fresh['pair_hashes'])):
                if not Path(str(path)).resolve().is_relative_to(upstream):
                    raise ValueError('Model points outside regenerated pair directory: '+str(path))
                check(f'regenerated pair integrity/{b}', digest(str(path)), str(expected))
            geometry = {key: fresh[key] for key in ('nside', 'zbin')}
        # Import the real production function, reusing the established catalog loader.
        spec = importlib.util.spec_from_file_location('analytic_noise', NOISE)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        params, ngal = module.catalog_noise_parameters(catalog, geometry)
        with np.load(args.reference_noise_ht) as old:
            check('catalog-derived noise variance', params['sigma_component_squared'], old['sigma_component_squared'], rtol=1e-12)
            check('selected galaxies', ngal, old['n_galaxies'])
        for label, path in [('signal', args.reference_signal_ht), ('noise', args.reference_noise_ht)]:
            with np.load(path) as old:
                check(label+'/reference model hash', str(old['model_sha256']), digest(reference_model))
                check(label+'/reference realizations', old['n_realizations'], 100)
                check(label+'/reference probability', old['keep_probability'], .001)
                check(label+'/reference base seed', old['base_seed'], 20260828)
        signal_out, noise_out = work/'signal_ht.npz', work/'noise_ht.npz'
        common = ['--model', model, '--mode', 'thinned', '--realizations', 100,
                  '--base-seed', 20260828, '--left-chunk', 0, '--right-chunk', 2000]
        run('signal-only HT', COMMON/'compute_covariance.py', *common, '--output', signal_out)
        # No --sims: calibration comes directly from the catalog, no realization arrays needed.
        run('signal plus shape-noise HT', NOISE, *common, '--catalog', catalog,
            '--signal-only-reference', signal_out, '--output', noise_out)
        keys = ('n_realizations','keep_probability','base_seed','theta','theta_edges',
                'full_counts','full_weight_sums','retained_counts','bin_seeds',
                'covariance_joint','covariance_sem','covariance_realizations')
        for label, newpath, oldpath in [('signal', signal_out, args.reference_signal_ht),
                                       ('noise', noise_out, args.reference_noise_ht)]:
            with np.load(newpath) as new, np.load(oldpath) as old:
                check(label+'/new model hash', str(new['model_sha256']), model_hash)
                for key in keys:
                    check(label+'/'+key, new[key], old[key], rtol=1e-12, atol=1e-20)
                if label == 'noise':
                    for key in ('pixel_noise_variance', 'sigma_component_squared', 'n_galaxies',
                                'covariance_blocks', 'covariance_noise_added'):
                        check(label+'/'+key, new[key], old[key], rtol=1e-12, atol=1e-20)
                    check('no simulation dependency', str(new['simulation_reference']), '')
                    check('no simulation count', new['simulation_nsims'], 0)
        state['outputs'] = {str(p.name): digest(p) for p in (signal_out, noise_out)}
        state['status'] = 'analytic_passed'
        state['all_comparisons_exact'] = all(c['exact'] for c in state['comparisons'])
        save()
        print('PASS: both analytic HT covariances and noise increment reproduced. '
              'Full simulation validation and comparison plots remain separate.', flush=True)
    except Exception as exc:
        state.update(status='failed', error=str(exc))
        save()
        raise


if __name__ == '__main__':
    main()
