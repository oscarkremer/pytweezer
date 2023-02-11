from pytweezer.shapes import Sphere
import numpy as np

radius = 1
shape = Sphere(radius)
#    [True, True, True, False]
rtp, n, ds = shape.boundary_points(npts=numr)
