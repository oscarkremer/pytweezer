import pytest
import numpy as np
from pytweezer.utils import rtp2xyz, xyz2rtp

def test_rtp2xyz():
    decimal = 3
    r = np.linspace(1, 11, 10)
    phi = np.pi/2*np.ones((10))
    theta = np.pi/2*np.ones((10))
    x, y, z = rtp2xyz(r, phi, theta)
    expected_x = np.zeros(10)
    expected_y = np.linspace(1, 11, 10)
    expected_z = np.zeros(10)
    np.testing.assert_array_almost_equal(x, expected_x, decimal=decimal, err_msg='Conversion from spherical to cartesian error in X')
    np.testing.assert_array_almost_equal(y, expected_y, decimal=decimal, err_msg='Conversion from spherical to cartesian error in Y')
    np.testing.assert_array_almost_equal(z, expected_z, decimal=decimal, err_msg='Conversion from spherical to cartesian error in Z')

def test_xyz2rtp():
    decimal = 3
    x = np.zeros(10)
    y = np.linspace(1, 11, 10)
    z = np.zeros(10)
    r, t, p = xyz2rtp(x, y, z)
    expected_r = np.linspace(1, 11, 10)
    expected_t = np.pi/2*np.ones(10)
    expected_p = np.pi/2*np.ones(10)
    np.testing.assert_array_almost_equal(r, expected_r, decimal=decimal, err_msg='Conversion from cartesian to spherical error in radius')
    np.testing.assert_array_almost_equal(t, expected_t, decimal=decimal, err_msg='Conversion from cartesian to spherical error in theta')
    np.testing.assert_array_almost_equal(p, expected_p, decimal=decimal, err_msg='Conversion from cartesian to spherical error in phi')