from pytweezer.shapes import Sphere
import numpy as np
from pytweezer.utils import angular_grid
from numpy import matlib as matlib

decimal = 3
radius = 1 
theta, phi = angular_grid(3, 3)
shape = Sphere(radius)
n = shape.normals(theta, phi)
sz = theta.size

print(n, matlib.repmat([1, 0, 0], sz, 1))
#    testCase.verifyThat(r, IsE