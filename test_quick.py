from pytweezer.utils import vswf_cart
import numpy as np



rtp = np.array([0,0,0])

A, B = vswf_cart(1, 0, rtp[0], rtp[1], rtp[2], 'regular')
print(A, B)
