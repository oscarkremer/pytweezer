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


def test_voxel_positions_with_offset():
    decimal = 3
    radius = 1 
    numr = 100
    sphere_1 = Sphere(radius)
    sphere_2 = Sphere(radius, position=np.array([1,0,0]))
    spacing = 0.1
    voxels1 = sphere_1.voxels(spacing, origin='world') + offset
    voxels2 = sphere_2.voxels(spacing, origin='world')


'''
  testCase.verifyThat(voxels2, IsEqualTo(voxels1, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Incorrect voxel positions');
'''