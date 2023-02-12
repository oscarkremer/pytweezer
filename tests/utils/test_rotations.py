import pytest
import numpy as np
from pytweezer.utils.rotations import rotx, roty, rotz


@pytest.mark.parametrize(
    ('input_matrix', 'expected'),
    (
        pytest.param(np.array([[1.0000, 0, 0],
                      [0, 0.7071, -0.7071],
                      [0, 0.7071, 0.7071]]), rotx(45), id='X rotation for 45 degrees'),
        pytest.param(np.array([[1, 0, 0],
                             [0, 0, -1], 
                             [0, 1, 0]]), rotx(90), id='X rotation for 90 degrees'),
        pytest.param(np.array([[1, 0, 0],
                              [0, -1, 0],
                              [0, 0, -1]]), rotx(180), id='X rotation for 180 degrees')
    )
)
def testRotx(input_matrix, expected):
  decimal = 3
  np.testing.assert_array_almost_equal(input_matrix, expected, decimal=decimal)


@pytest.mark.parametrize(
    ('input_matrix', 'expected'),
    (
        pytest.param(np.array([[0.7071, 0, 0.7071],
                      [0, 1.000, 0.],
                      [-0.7071, 0, 0.7071]]), roty(45), id='Y rotation for 45 degrees'),
        pytest.param(np.array([[0, 0, 1],
                             [0, 1, 0], 
                             [-1, 0, 0]]), roty(90), id='Y rotation for 90 degrees'),
        pytest.param(np.array([[-1, 0, 0],
                              [0, 1, 0],
                              [0, 0, -1]]), roty(180), id='Y rotation for 180 degrees')
    )
)
def testRoty(input_matrix, expected):
  decimal = 3
  np.testing.assert_array_almost_equal(input_matrix, expected, decimal=decimal)


@pytest.mark.parametrize(
    ('input_matrix', 'expected'),
    (
        pytest.param(np.array([[0.7071, -0.7071, 0],
                      [0.7071, 0.7071, 0],
                      [0, 0, 1.000]]), rotz(45), id='Z rotation for 45 degrees'),
        pytest.param(np.array([[0, -1, 0],
                             [1, 0, 0], 
                             [0, 0, 1]]), rotz(90), id='Z rotation for 90 degrees'),
        pytest.param(np.array([[-1, 0, 0],
                              [0, -1, 0],
                              [0, 0, 1]]), rotz(180), id='Z rotation for 180 degrees')
    )
)
def testRotz(input_matrix, expected):
  decimal = 3
  np.testing.assert_array_almost_equal(input_matrix, expected, decimal=decimal)