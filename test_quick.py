from pytweezer.shapes import Sphere
import numpy as np

radius = 1
shape = Sphere(radius)
numr=100
#    [True, True, True, False]
r, t, p, n_rho, n_theta, n_phi, ds = shape.boundary_points(npts=numr)
