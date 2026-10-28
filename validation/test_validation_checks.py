"""Directly runnable checks for result comparison; no pytest or astronomy packages."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np

from run_from_catalog import compare_array, matched_rows

path = Path(__file__).resolve().parents[1]/'examples/xmm_nside512/shear_XMM_Nside512_signal_plus_shape_noise_sims/simulate_shear_shape_noise.py'
spec = importlib.util.spec_from_file_location('noise', path)
noise = importlib.util.module_from_spec(spec)
spec.loader.exec_module(noise)


class ValidationChecks(unittest.TestCase):
    def test_counts_and_uint64_seed_changes_are_not_rounded_away(self):
        a = np.array([2**63], dtype=np.uint64)
        b = np.array([2**63+1], dtype=np.uint64)
        self.assertFalse(compare_array(a, b)['passed'])
        self.assertFalse(compare_array([[0, 1]], [[1, 0]])['passed'])

    def test_nonfinite_or_wrong_shape_is_never_passed(self):
        for a, b in [([np.nan], [np.nan]), ([np.inf], [np.inf]), ([[1]], [1])]:
            self.assertFalse(compare_array(a, b)['passed'])

    def test_small_covariance_differences_are_not_hidden_by_numpy_default_atol(self):
        self.assertFalse(compare_array([1e-12], [2e-12], rtol=2e-9, atol=1e-20)['passed'])
        result = compare_array([1e-12+1e-22], [1e-12], rtol=2e-9, atol=1e-20)
        self.assertTrue(result['passed'])
        self.assertFalse(result['exact'])

    def test_seed_alignment_and_missing_or_duplicate_seed(self):
        self.assertEqual(matched_rows([2, 0], [0, 1, 2]), [2, 0])
        with self.assertRaises(KeyError):
            matched_rows([3], [0, 1, 2])
        with self.assertRaises(ValueError):
            matched_rows([0], [0, 0])

    def test_catalog_noise_calibration_and_independent_rng(self):
        e1, e2, w = np.array([.1, -.3, .2]), np.array([.2, .1, -.4]), np.array([1., 2., 3.])
        params = noise.noise_parameters(e1, e2, w)
        m1, m2 = np.average(e1**2, weights=w), np.average(e2**2, weights=w)
        expected = (m1/(4*(1-m1)**2) + m2/(4*(1-m2)**2))/2
        np.testing.assert_allclose(params['sigma_component_squared'], expected, rtol=1e-15)
        np.random.seed(50)
        before = np.random.get_state()
        generated = noise.draw_noise(100000, expected, 13)
        after = np.random.get_state()
        np.testing.assert_array_equal(before[1], after[1])
        self.assertEqual(before[2:], after[2:])
        np.testing.assert_array_equal(generated, noise.draw_noise(100000, expected, 13))
        np.testing.assert_allclose(generated.var(axis=1), expected, rtol=.02)
        self.assertLess(abs(np.corrcoef(generated)[0, 1]), .02)


if __name__ == '__main__':
    unittest.main()
