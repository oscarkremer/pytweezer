import pytest
import warnings
import numpy as np
warnings.filterwarnings("ignore", category=PendingDeprecationWarning)
from numpy import matlib as matlib
from pytweezer.shapes import Sphere
from pytweezer.utils import angular_grid

def test_sphere_radii():
    decimal = 3
    radius = 1 
    theta, phi = angular_grid(3, 3)
    shape = Sphere(radius)
    r = shape.radii(theta, phi)
    sz = theta.size
    np.testing.assert_array_almost_equal(r, matlib.repmat(1, 1, sz).reshape((sz)), decimal=decimal,
        err_msg='Error for radii function of Sphere shape')

def test_sphere_normals():
    decimal = 3
    radius = 1 
    theta, phi = angular_grid(3, 3)
    shape = Sphere(radius)
    n = shape.normals(theta, phi)
    sz = theta.size
    np.testing.assert_array_almost_equal(n, matlib.repmat([1, 0, 0], sz, 1), decimal=decimal,
        err_msg='Error for normal function of Sphere shape')

def test_sphere_axial_symmetry():
    decimal = 3
    radius = 1 
    theta, phi = angular_grid(3, 3)
    shape = Sphere(radius)
    sz = theta.size
    _, _, rotsym = shape.axial_symmetry()
    assert rotsym==0, 'Error for axial symmetry of Sphere shape'