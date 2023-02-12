import pytest
import numpy as np
from pytweezer.shapes import Sphere

def test_inside_xyz_boundary_points_values_r():  
    radius = 1.0
    offset = np.array([0, 0, 2])
    shape = Sphere(radius, position=offset)
    z = np.linspace(-10, 10, 100).reshape((100,1))
    mask = shape.inside_xyz(0, 0, z, origin='world') 
    assert mask.shape==z.shape, 'Error in dimension of returned mask'


'''
def test_voxel_positions_with_offset():
    decimal = 3
    radius = 1 
    numr = 100
    shape = Sphere(radius)
    _, _, _, n_rho, _, _, _ = shape.boundary_points(npts=numr)
    np.testing.assert_array_almost_equal(n_rho, np.ones((n_rho.shape)), 
        decimal=decimal, err_msg='Error computing n_rho vector in boundary_points')
          import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-6;

  offset = [1;0;0];
  radius = 1.0;
  shape1 = ott.shapes.Sphere(radius);
  shape2 = ott.shapes.Sphere(radius, offset);

  spacing = 0.1;
  voxels1 = shape1.voxels(spacing, 'origin', 'world') + offset;
  voxels2 = shape2.voxels(spacing, 'origin', 'world');

  testCase.verifyThat(voxels2, IsEqualTo(voxels1, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect voxel positions');
'''