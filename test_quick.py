from pytweezer.shapes import Sphere
import numpy as np
from pytweezer.beams import Gaussian
from pytweezer.utils import angular_grid, vswf_cart
from numpy import matlib as matlib


beam = Gaussian()
beam.power = 1.0
dz = np.pi/2
tol = 1e-6
beam.translate_z(z=dz)
