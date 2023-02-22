from pytweezer.beams import Gaussian
from pytweezer.t_matrix import TMatrixMie

beam = Gaussian(power=1.0)

T = TMatrixMie(radius=1.0, index_r=1.2)


