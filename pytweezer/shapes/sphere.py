class Sphere(StarShape, AxisymShape):

    def __init__(self, radius, position=[]):
        self.radius = radius
        if position:
            self.position = position

    def get_maxRadius(self, shape):
        return shape.radius

    def get_volume(self, shape):
      v = 4/3*np.pi*np.power(shape.radius, 3)

    def get_perimiter(self, shape):
      return 2.0 * np.pi * shape.radius

    def radii(self, shape, theta, phi):
        theta, phi = matchsize(theta, phi)
        return np.ones(theta.shape)*shape.radius

    def normals(self, shape, theta, phi):
        theta, phi = matchsize(theta, phi)
        return np.ones(theta.shape) * np.array([ 1, 0, 0 ])

    def boundarypoints(shape, **kwargs):
        ntheta = shape.boundarypoints_npts(**kwargs)
        theta = np.arange(0.0, (np.pi/(ntheta-1)):np.pi)
        phi = np.zeros(theta.shape)
        xyz = shape.locations(theta, phi)
        nxyz = xyz/shape.radius
        n, rtp = xyzv2rtpv(nxyz, xyz)
        ds = boundarypoints_area(xyz[:, 0], xyz[:, 2], xyz[:, 0], xyz[:, 2], rtp)
        return rtp, n, ds

    def axial_symmetry(self, shape):
        return 0, 0, 0