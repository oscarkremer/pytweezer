import pytest
import numpy as np
from pytweezer.utils import paraxial_beam_waist


@pytest.mark.parametrize(
    ('input_x', 'expected'),
    (
        pytest.param(0, 1, id='Paraxial beam waist for input = 0'),
        pytest.param(1, 1.5009, id='Paraxial beam waist for input = 1'),
        pytest.param(2, 1.7739, id='Paraxial beam waist for input = 2')
    )
)
def test_paraxial_beam_waist(input_x, expected):
    decimal = 3
    w = paraxial_beam_waist(input_x)
    np.testing.assert_array_almost_equal(w, expected, decimal=decimal, err_msg='Waist failed test for equality')
