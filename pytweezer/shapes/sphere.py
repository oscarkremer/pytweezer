import numpy as np
from .axi_sym_shape import AxiSymShape
from .star_shape import StarShape
from pytweezer.utils import angular_grid, match_size, xyz2rtp, xyzv2rtp, rtp2xyz

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
        ntheta = self.boundary_points_npts(**kwargs)
        theta = np.arange(0.0, np.pi+np.pi/(ntheta-1), np.pi/(ntheta-1))
        phi = np.zeros(theta.shape)
        x, y, z = self.locations(theta, phi)
        nx, ny, nz = x/self.__radius__, y/self.__radius__, z/self.__radius__
        n, r, t, p = xyzv2rtpv(nx, ny, nz, x, y, z)
        ds = boundary_points_area(x, y, z, z, rtp)
        return rtp, n, ds

    def axial_symmetry(self, shape):
        return 0, 0, 0