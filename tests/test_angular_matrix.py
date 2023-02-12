import pytest
import numpy as np
from pytweezer.utils import angular_grid

@pytest.mark.parametrize(
    ('input_x', 'expected', 'style'),
    (
        pytest.param((3, 3), (np.array([0.52359878, 1.57079633, 2.61799388,
                                       0.52359878, 1.57079633, 2.61799388,
                                       0.52359878, 1.57079633,2.61799388]),
                            np.array([0, 0, 0, 2.0943951, 2.0943951, 2.0943951,
                                    4.1887902, 4.1887902, 4.188790])), 'column',  id='Column case for squared matrix'),
        pytest.param((3, 3), (np.array([[0.5236, 0.5236, 0.5236],
                                        [1.5708, 1.5708, 1.5708],
                                        [2.6180, 2.6180, 2.6180]]),
                            np.array([[0, 2.0944, 4.1888],
                                     [0, 2.0944, 4.1888],
                                     [0, 2.0944, 4.1888]])), 'matrix',  id='Matrix case for squared matrix'),
        pytest.param((3, 3), (np.array([[0.5236], [1.5708], [2.6180],
                                        [0.5236], [1.5708], [2.6180],
                                        [0.5236], [1.5708], [2.6180]]),
                            np.array([[0], [0], [0], 
                                      [2.0944], [2.0944], [2.0944], 
                                      [4.1888], [4.1888], [4.1888]])), 'points',  id='Points case for squared matrix'),
        pytest.param((2, 3), (np.array([0.7854, 2.3562, 0.7854, 2.3562, 0.7854, 2.3562]),
                            np.array([0, 0, 2.0944, 2.0944, 4.1888, 4.1888])), 'column',  id='Column case for non-squared matrix'),
        pytest.param((2, 3), (np.array([[0.7854, 0.7854, 0.7854],
                                        [2.3562, 2.3562, 2.3562]]),
                            np.array([[0, 2.0944, 4.1888],
                                     [0, 2.0944, 4.1888]])), 'matrix',  id='Matrix case for non-squared matrix'),
        pytest.param((2, 3), (np.array([[0.7854], [2.3562],
                                        [0.7854], [2.3562],
                                        [0.7854], [2.3562]]),
                            np.array([[0], [0], 
                                      [2.0944], [2.0944], 
                                      [4.1888], [4.1888]])), 'points',  id='Points case for non-squared matrix'),
        pytest.param((2, 3), (np.array([0.5236, 1.5708, 2.6180]),
                            np.array([0, 2.0944, 4.1888])), 'row',  id='Value error case for using a style which is not allowed'),
    )
)
def test_angular_matrix(input_x, expected, style):
    decimal = 3
    if style not in ('column', 'matrix', 'points'):
        with pytest.raises(ValueError):
            theta, phi = angular_grid(input_x[0], input_x[1], style)
    else:
        theta, phi = angular_grid(input_x[0], input_x[1], style)
        np.testing.assert_array_almost_equal(theta, expected[0], decimal=decimal, err_msg='Theta failed test for equality')
        np.testing.assert_array_almost_equal(phi, expected[1], decimal=decimal, err_msg='Phi failed test for equality')