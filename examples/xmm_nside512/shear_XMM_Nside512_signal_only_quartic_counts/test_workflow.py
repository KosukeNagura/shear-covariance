"""Small synthetic integration test for the thinned covariance path."""
import os
from pathlib import Path
import subprocess
import sys

import healpy as hp
import numpy as np

import shear_core as core


ROOT = Path(__file__).resolve().parent


def test_thinned_cli_on_small_synthetic_model(tmp_path):
    pixels = np.array([110, 207, 330, 404])
    vectors = np.array(hp.pix2vec(8, pixels)).T
    weights = np.array([1.0, 2.0, 0.7, 1.3])
    pairs_by_bin = [np.array([[0, 1], [0, 2]]), np.array([[1, 2], [2, 3]])]
    pair_paths = []
    pair_hashes = []
    weight_sums = []
    for b, pairs in enumerate(pairs_by_bin):
        path = tmp_path / f'pairs_bin_{b:02d}.npy'
        np.save(path, pairs)
        pair_paths.append(str(path))
        pair_hashes.append(core.file_hash(path))
        weight_sums.append(np.sum(weights[pairs[:, 0]] * weights[pairs[:, 1]]))

    cl = np.array([0.0, 0.0, 0.7, 0.4, 0.2])
    theta_grid = np.linspace(0, np.pi, 60001)
    xip, xim = core.xi_from_cl(theta_grid, cl)
    signal = core.signal_covariance(vectors, theta_grid, xip, xim, xip[0])
    model = tmp_path / 'model.npz'
    core.exclusive_save(
        model, vectors=vectors, weights=weights, pix_idx=pixels, nside=8, zbin=4,
        cl_ee=cl, signal=signal, pair_paths=np.asarray(pair_paths),
        pair_hashes=np.asarray(pair_hashes), full_weight_sums=np.asarray(weight_sums),
        full_counts=np.array([2, 2]), theta=np.array([10.0, 50.0]),
        theta_edges=np.array([1.0, 30.0, 179.0]), estimator='synthetic test',
    )

    output = tmp_path / 'ht_p0.001.npz'
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    result = subprocess.run(
        [sys.executable, str(ROOT / 'compute_covariance.py'), '--model', str(model),
         '--mode', 'thinned', '--realizations', '4', '--output', str(output)],
        env=env, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    with np.load(output) as data:
        assert int(data['n_realizations']) == 4
        assert float(data['keep_probability']) == 0.001
        assert data['covariance_joint'].shape == (4, 4)
        assert data['covariance_sem'].shape == (4, 4)
        assert np.all(np.isfinite(data['covariance_joint']))
        assert np.all(np.isfinite(data['covariance_sem']))

    # The command must preserve the existing no-overwrite safeguard.
    rerun = subprocess.run(
        [sys.executable, str(ROOT / 'compute_covariance.py'), '--model', str(model),
         '--mode', 'thinned', '--realizations', '4', '--output', str(output)],
        env=env, capture_output=True, text=True,
    )
    assert rerun.returncode != 0
    assert 'Refusing existing final output' in rerun.stderr
