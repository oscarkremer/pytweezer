import pytest
import numpy as np
from pytweezer.utils.three_wide import three_wide


@pytest.mark.parametrize(
    ('expected', 'input_matrix'),
    (
        pytest.param(np.array([[1, 1, 1]]), three_wide(1), id='Three Wide with integer input'),
        pytest.param(np.array([[1.0, 1.0, 1.0]]), three_wide(1.0), id='Three Wide with float input'),
        pytest.param(np.array([[1, 1, 1],
                              [2, 2, 2]]), three_wide(np.array([1, 2])), id='Three Wide with array input')
    )
)
def test_three_wide_computation(input_matrix, expected):
  decimal = 3
  np.testing.assert_array_almost_equal(input_matrix, expected, decimal=decimal)


def test_three_wide_argument_type():
    with pytest.raises(TypeError):
        three_wide('test')

def test_three_wide_argument_value():
    with pytest.raises(ValueError):
        three_wide(np.array([[1, 1, 1]]))