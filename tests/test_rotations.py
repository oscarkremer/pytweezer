import pytest
import numpy as np
from pytweezer.utils.rotations import rotx, roty, rotz
#@pytest.fixture
#def cart():
#    return ShoppingCart(5)
#

#def square(n: int) -> int:
#    return n*n

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

'''
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
        square(input_x)v
'''

'''
def testRotz(testCase)

  tol = 1.0e-4;

  rotz45 = [0.7071   -0.7071         0
            0.7071    0.7071         0
                 0         0    1.0000];
               
  rotz90= [0, -1, 0; 1, 0, 0; 0, 0, 1];
  rotz180= [-1, 0, 0; 0, -1, 0; 0, 0, 1];

  val = ott.utils.rotz(45);
  testCase.verifyThat(val, IsEqualTo(rotz45, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect rotz values');

  val = ott.utils.rotz([1, 2, 3; 4, 5, 6]);
  testCase.verifyThat(size(val), IsEqualTo([3, 3*6]), ...
      'Incorrect size for rotz matrix input');

  val = ott.utils.rotz([1, 2, 3; 4, 5, 6], 'usecell', true);
  testCase.verifyThat(size(val), IsEqualTo([2, 3]), ...
      'Incorrect size for rotz matrix input cell output');

  val = ott.utils.rotz([45, 90; 0, 180]);
  res = [rotz45, eye(3), rotz90, rotz180];
  testCase.verifyThat(val, IsEqualTo(res, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect output for rotz matrix input');

  val = ott.utils.rotz([45, 90; 0, 180], 'usecell', true);
  res = {rotz45, rotz90; eye(3), rotz180};
  testCase.verifyThat(val, IsEqualTo(res, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect output for rotz cell input');


def testRoty(testCase)
  tol = 1.0e-4;

  roty45 = [0.7071         0    0.7071
                 0    1.0000         0
           -0.7071         0    0.7071];
               
  roty90= [0, 0, 1; 0, 1, 0; -1, 0, 0];
  roty180= [-1, 0, 0; 0, 1, 0; 0, 0, -1];

  val = ott.utils.roty(45);
  testCase.verifyThat(val, IsEqualTo(roty45, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect roty values');

  val = ott.utils.roty([1, 2, 3; 4, 5, 6]);
  testCase.verifyThat(size(val), IsEqualTo([3, 3*6]), ...
      'Incorrect size for roty matrix input');

  val = ott.utils.roty([1, 2, 3; 4, 5, 6], 'usecell', true);
  testCase.verifyThat(size(val), IsEqualTo([2, 3]), ...
      'Incorrect size for roty matrix input cell output');

  val = ott.utils.roty([45, 90; 0, 180]);
  res = [roty45, eye(3), roty90, roty180];
  testCase.verifyThat(val, IsEqualTo(res, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect output for roty matrix input');

  val = ott.utils.roty([45, 90; 0, 180], 'usecell', true);
  res = {roty45, roty90; eye(3), roty180};
  testCase.verifyThat(val, IsEqualTo(res, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect output for roty cell input');
'''
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