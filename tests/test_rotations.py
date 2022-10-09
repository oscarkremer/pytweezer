import pytest

#@pytest.fixture
#def cart():
#    return ShoppingCart(5)
#

def square(n: int) -> int:
    return n*n

#@pytest.mark.parametrize(
#    ('input_x', 'expected'),
#    (
#        pytest.param(1, 1, id='trivial positive case'),
#        pytest.param(-1, 1, id='trivial negative case'),
#        pytest.param(2, 4, id='positive case')
#    )
#)
#def test_square(input_x, expected, cart):
#    assert square(input_x) == expected


@pytest.mark.parametrize(
    'input_x',
    (
        'a',
        [],
        (),
    )
)
def test_square_error(input_x):
    with pytest.raises(TypeError):
        square(input_x)