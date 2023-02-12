import pytest
import numpy as np
from pytweezer.utils import sbesselj, sbesselh

@pytest.mark.parametrize(
    ('n', 'kr', 'verbose', 'expected_jn', 'expected_djn'),
    (
        pytest.param(1, 1e-16, False, np.array([[3.3333e-17]]), np.array([[-5.5511e-17]]), 
            id='Spherical Bessel Function with n and k being numbers'),
        pytest.param(np.array([1, 2, 3]), 1, True, np.array([[0.3012, 0.0620, 0.0090]]), 
        np.array([[0.5403, 0.1771, 0.0350]]), id='Spherical Bessel Function with n being an array.'),

        pytest.param(np.array([1, 2, 3]), np.array([np.pi/10, np.pi/2, np.pi]),  False,
                np.array([[0.1037, 0.0065, 0.0003],
                [0.4053, 0.1374, 0.0321],
                [0.3183, 0.3040, 0.1655]]), 
                np.array([[0.6536, 0.0621, 0.0037],
                [0.3786, 0.2303, 0.0761],
                [-0.1013, 0.1248, 0.1460]]), id='Spherical Bessel Function with n and k being an array.'),
    )
)
def test_spherical_bessel_function(n, kr, verbose, expected_jn, expected_djn):
    decimal = 3
    jn, djn = sbesselj(n, kr, verbose=verbose)
    np.testing.assert_array_almost_equal(jn, expected_jn, decimal=decimal, 
        err_msg='Error for spherical bessel function computation!')
    np.testing.assert_array_almost_equal(djn, expected_djn, decimal=decimal,
        err_msg='Error for spherical bessel function derivative computation')


@pytest.mark.parametrize(
    ('n', 'kr'),
    (
        pytest.param('test-input', 1e-16, 
            id='Spherical Bessel Function with n begin string'),
        pytest.param('test-input', 1e-16, id='Spherical Bessel Function with n being list.'),
    )
)
def test_spherical_bessel_args(n, kr):
    with pytest.raises(TypeError):
        sbesselj(n, kr)