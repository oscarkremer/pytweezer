import numpy as np
from pytweezer.utils import vswf_cart
from pytweezer.utils.rotations import rotx, rotz

def test_vswf_cart():
    decimal = 3
    A, B = vswf_cart(1, 0, 0, 0, 0, htype='regular')
    assert np.isfinite(A).all() and np.isfinite(B).all(), 'A and B contain only finite elements'