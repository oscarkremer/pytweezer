import pytest
import numpy as np
from pytweezer.utils import sbesselj, sbesselh

@pytest.mark.parametrize(
    ('n', 'kr', 'expected_jn', 'expected_djn'),
    (
        pytest.param(np.array([1, 2, 3]), 1, np.array([[0.3012, 0.0620, 0.0090]]), 
        np.array([[0.5403, 0.1771, 0.0350]]), id='Spherical Bessel Function with n being an array.'),
    )
)
def test_spherical_bessel_function(n, kr, expected_jn, expected_djn):
  decimal = 3
  jn, djn = sbesselj(n, kr)
  np.testing.assert_array_almost_equal(jn, expected_jn, decimal=decimal, 
    err_msg='Error for spherical bessel function computation!')
  np.testing.assert_array_almost_equal(djn, expected_djn, decimal=decimal,
    err_msg='Error for spherical bessel function derivative computation')

