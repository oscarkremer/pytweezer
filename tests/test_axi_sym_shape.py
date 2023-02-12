import pytest
import numpy as np
from pytweezer.shapes import Sphere

def test_sphere_boundary_points():
    decimal = 3
    radius = 1 
    shape = Sphere(radius)
    np.testing.assert_array_almost_equal(shape.perimeter, 6.2832, decimal=decimal, err_msg='Error computing perimeter of sphere')
    numr = 100
    r, t, p, n_rho, n_theta, n_phi, ds = shape.boundary_points(npts=numr)
'''
  testCase.verifyThat(length(rtp(:, 1)), IsEqualTo(numr), ...
      'Incorrect number of values returned');

  testCase.verifyThat(rtp(:, 1), IsEqualTo(ones(size(rtp(:, 1))), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect values for r');

  testCase.verifyThat(n(:, 1), IsEqualTo(ones(size(n(:, 1))), ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect values for r');
def test_sphere_volume():
    decimal = 3
    radius = 1 
    shape = Sphere(radius)
    np.testing.assert_array_almost_equal(shape.volume, 4.1888, decimal=decimal, err_msg='Error computing volume of sphere')


def test_sphere_inside():
    radius = 1
    shape = Sphere(radius)
    x = np.append(0.5*1*np.random.rand(1, 3), 4.0)
    y = np.append(0.5*1*np.random.rand(1, 3), 4.0)
    z = np.append(0.5*1*np.random.rand(1, 3), 4.0)
    for value, expected in zip(shape.inside_xyz(x,y,z), [True, True, True, False]):
        assert value==expected, 'Error in inside_xyz method for Spherical Shape'
'''    