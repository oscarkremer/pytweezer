from pytweezer.shapes import Sphere
import numpy as np
from pytweezer.beams import Gaussian
from pytweezer.beams.translation import translate_z, translate
from pytweezer.utils import angular_grid, vswf_cart
from numpy import matlib as matlib
from copy import copy

beam = Gaussian()

#dz=np.pi/2
#beam_2, A, B = translate_z(copy(beam), z=dz)
#print(beam_2.a)