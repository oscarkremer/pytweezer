import numpy as np
from numba import njit
from numba.pycc import CC
cc = CC('meshgrid')

@njit(nopython=True, cache=True)
@cc.export('meshgrid', 'Tuple((f8[:,:], i8[:,:]))(f8[:], i8[:])')
def meshgrid(x, y):
    xx = np.empty(shape=(y.size, x.size), dtype=x.dtype)
    yy = np.empty(shape=(y.size, x.size), dtype=y.dtype)
    for j in range(y.size):
        for k in range(x.size):
            xx[j, k] = x[k]
            yy[j, k] = y[j]
    return xx, yy

if __name__ == '__main__':
    cc.compile()
