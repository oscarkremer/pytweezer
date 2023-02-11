import pytest
import numpy as np
from pytweezer.shapes import Shape

@pytest.mark.parametrize(
    ('point', 'radius', 'expected'),
    (
          x = [0.5*0.5*rand(1, 3), 4.0];

          id='Inside point case 1'),
    )
)
def test_sphere_shape(point, radius, expected):
    decimal = 3
    shape = Shape.simple('sphere', r=radius)
    if style not in ('column', 'matrix', 'points'):
        with pytest.raises(ValueError):
            theta, phi = angular_grid(input_x[0], input_x[1], style)
    else:
        theta, phi = angular_grid(input_x[0], input_x[1], style)
        np.testing.assert_array_almost_equal(theta, expected[0], decimal=decimal, err_msg='Theta failed test for equality')
        np.testing.assert_array_almost_equal(phi, expected[1], decimal=decimal, err_msg='Phi failed test for equality')