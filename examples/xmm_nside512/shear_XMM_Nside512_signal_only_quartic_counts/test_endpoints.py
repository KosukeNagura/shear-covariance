"""Endpoint and low-multipole conventions for the spin-2 correlation."""
import numpy as np
import pytest

from shear_core import xi_from_cl


def test_zero_lag_and_ell_zero_one_do_not_affect_spin2_correlation():
    cl = np.array([9.0, 8.0, 0.3, 0.2, 0.1])
    theta = np.array([0.0, 0.05, 0.2, np.pi])
    xip, xim = xi_from_cl(theta, cl)

    expected_zero = np.sum((2*np.arange(2, len(cl)) + 1) * cl[2:]) / (4*np.pi)
    np.testing.assert_allclose(xip[0], expected_zero, rtol=1e-15)
    assert xim[0] == 0.0

    cleaned = cl.copy()
    cleaned[:2] = 0.0
    clean_xip, clean_xim = xi_from_cl(theta, cleaned)
    np.testing.assert_array_equal(xip, clean_xip)
    np.testing.assert_array_equal(xim, clean_xim)


def test_bad_theta_is_rejected():
    with pytest.raises(ValueError, match=r'\[0, pi\]'):
        xi_from_cl(np.array([-0.01]), np.array([0.0, 0.0, 1.0]))
