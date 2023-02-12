import pytest
import numpy as np
from pytweezer.shapes import Sphere

def test_sphere_boundary_points():
    decimal = 3
    radius = 1 
    numr = 100
    shape = Sphere(radius)
    r, t, p, n_rho, n_theta, n_phi, ds = shape.boundary_points(npts=numr)
    np.testing.assert_array_almost_equal(shape.perimeter, 6.2832, decimal=decimal, err_msg='Error computing perimeter of sphere')

def test_sphere_boundary_points_shape_rtp():
    decimal = 3
    radius = 1 
    numr = 100
    shape = Sphere(radius)
    r, _, _, _, _, _, _ = shape.boundary_points(npts=numr)
    assert r.shape[0]==numr, 'Error in the shape of the returned vector for radial coordinate'

def test_sphere_boundary_points_values_r():
    decimal = 3
    radius = 1 
    numr = 100
    shape = Sphere(radius)
    r, _, _, _, _, _, _ = shape.gboundary_points(npts=numr)
    np.testing.assert_array_almost_equal(r, np.ones((r.shape)), 
        decimal=decimal, err_msg='Error computing r vector in boundary_points')

def test_sphere_boundary_points_values_n():
    decimal = 3
    radius = 1
    numr = 100
    shape = Sphere(radius)
    _, _, _, n_rho, _, _, _ = shape.boundary_points(npts=numr)
    np.testing.assert_array_almost_equal(n_rho, np.ones((n_rho.shape)), 
        decimal=decimal, err_msg='Error computing n_rho vector in boundary_points')
