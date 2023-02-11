import numpy as np
from .axi_sym_shape import AxiSymShape
from .star_shape import StarShape
from pytweezer.utils import angular_grid, match_size, xyz2rtp, rtp2xyz

class Sphere(StarShape, AxiSymShape):

    def __init__(self, radius, position=[]):
        self.__radius__ = radius
        if position:
            self.position = position
        self.perimeter = self.get_perimeter()
        self.volume = self.get_volume()

    def get_max_radius(self):
        return self.__radius__

    def get_volume(self):
        return 4/3*np.pi*np.power(self.__radius__, 3)

    def get_perimeter(self):
      return 2.0 * np.pi * self.__radius__

    def radii(self, theta, phi):
        theta, phi = match_size(theta, phi)
        print(theta)
        return np.ones(theta.shape)*self.__radius__

    def normals(self, theta, phi):
        theta, phi = match_size(theta, phi)
        return np.ones(theta.shape) * np.array([ 1, 0, 0 ])

    def boundary_points(self, **kwargs):
        ntheta = self.boundary_points_npts(kwargs)
        theta = np.arange(0.0, (np.pi/(ntheta-1)),np.pi)
        phi = np.zeros(theta.shape)
        xyz = shape.locations(theta, phi)
        nxyz = xyz/self.__radius__
        n, rtp = xyzv2rtpv(nxyz, xyz)
        ds = boundary_points_area(xyz[:, 0], xyz[:, 2], xyz[:, 0], xyz[:, 2], rtp)
        return rtp, n, ds

    def axial_symmetry(self, shape):
        return 0, 0, 0