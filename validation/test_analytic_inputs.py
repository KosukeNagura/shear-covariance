"""Small numeric input-routing tests; stdlib/NumPy only, no production data."""
import contextlib
import importlib.util
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from run_analytic_validation import upstream_ready

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT/'examples/xmm_nside512/shear_XMM_Nside512_signal_plus_shape_noise_quartic_counts/compute_covariance_noise.py'
spec = importlib.util.spec_from_file_location('analytic', PATH)
analytic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analytic)


class AnalyticInputs(unittest.TestCase):
    def test_subtraction_checks_model_and_sampling_and_preserves_value(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'signal.npz'
            seeds = np.array([[2**63+1]], dtype=np.uint64)
            counts = np.array([[3]])
            base = np.array([[1., .2], [.2, 2.]])*1e-11
            total = base*1.5
            np.savez(path, covariance_joint=base, bin_seeds=seeds, retained_counts=counts,
                     keep_probability=.001, base_seed=20260828, model_sha256='same-model')
            args = (path, total, seeds, counts, .001, 20260828)
            np.testing.assert_array_equal(analytic.checked_noise_increment(*args, 'same-model'), total-base)
            with self.assertRaises(ValueError):
                analytic.checked_noise_increment(*args, 'different-model')
            with self.assertRaises(ValueError):
                analytic.checked_noise_increment(path, total, seeds+np.uint64(1), counts, .001, 20260828)
            with self.assertRaises(ValueError):
                analytic.checked_noise_increment(path, total, seeds, counts+1, .001, 20260828)

    def test_catalog_route_and_archive_route_produce_same_covariance(self):
        # Mock only FITS/HEALPix I/O boundaries. The actual covariance kernels run.
        with tempfile.TemporaryDirectory() as temp:
            d = Path(temp)
            rng = np.random.default_rng(71)
            vectors = rng.normal(size=(4,3))
            vectors /= np.linalg.norm(vectors, axis=1)[:,None]
            weights = np.array([1.,2.,.7,1.3])
            pairs = np.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
            paths = []
            for b in range(10):
                p = d/f'pairs_{b}.npy'; np.save(p, pairs); paths.append(str(p))
            a = rng.normal(size=(8,8))*1e-3
            signal = (a@a.T).reshape(4,2,4,2).transpose(0,2,1,3)
            theta = np.arange(1,11,dtype=float)
            model = d/'model.npz'
            np.savez(model, signal=signal, vectors=vectors, weights=weights,
                     pix_idx=np.arange(4), nside=8, zbin=4, theta=theta,
                     theta_edges=np.arange(11,dtype=float), full_counts=np.full(10,6),
                     full_weight_sums=np.full(10,np.sum(weights[pairs[:,0]]*weights[pairs[:,1]])),
                     pair_paths=paths, pair_hashes=[analytic.file_hash(p) for p in paths])
            sims = d/'sims.npz'
            np.savez(sims, sigma_component_squared=.067, n_galaxies=9, nsims=3,
                     nside=8, zbin=4, theta=theta, covariance_total=np.eye(20),
                     covariance_signal=np.eye(20)*2, covariance_noise=np.eye(20)*3)
            outputs = []
            for label, extras in [('catalog', []), ('archive', ['--sims', str(sims)])]:
                out = d/(label+'.npz')
                # Full mode only on this four-pixel fixture avoids an empty HT sample.
                args = [str(PATH), '--model',str(model),'--catalog',str(d/'unused.fits'),
                        '--mode','full','--output',str(out), *extras]
                with patch.object(sys, 'argv', args), contextlib.redirect_stdout(io.StringIO()), \
                     patch.object(analytic, 'catalog_noise_parameters', return_value=({'sigma_component_squared':.067},9)), \
                     patch.object(analytic, 'read_pixel_noise', return_value=np.array([.01,.02,.03,.04])):
                    analytic.main()
                outputs.append(out)
            with np.load(outputs[0]) as cat, np.load(outputs[1]) as sim:
                for key in ('covariance_joint','covariance_sem','covariance_blocks','covariance_realizations',
                            'pixel_noise_variance','sigma_component_squared','n_galaxies','retained_counts','bin_seeds'):
                    np.testing.assert_array_equal(cat[key],sim[key],err_msg=key)
                self.assertEqual(str(cat['noise_calibration_source']), 'catalog')
                self.assertEqual(str(cat['simulation_reference']), '')
                self.assertEqual(cat['simulation_covariance_total'].size, 0)
                self.assertGreater(np.max(np.abs(cat['covariance_joint'])), 0.)

    def test_partial_upstream_report_cannot_pass(self):
        with tempfile.TemporaryDirectory() as temp:
            d=Path(temp)
            (d/'validation_report.json').write_text('{"comparisons": []}')
            with self.assertRaises(ValueError):
                upstream_ready(d)


if __name__ == '__main__':
    unittest.main()
