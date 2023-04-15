import numpy as np
from scipy.linalg import block_diag
from numba import jit

#@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
def wigner_rotation_matrix(nmax, R):
    '''    % WIGNER_ROTATION_MATRIX rotation matrix for rotation of spherical
    % harmonics or T-matrices.
    %
    % D = WIGNER_ROTATION_MATRIX(nmax,R) calculates the rotation matrix
    % for the VSH given a 3x3 coordinate rotation matrix R.  Usage: a' = D a.
    %
    % This method from Choi et al., J. Chem. Phys. 111: 8825-8831 (1999)
    % Note change in notation - here, use x' = Rx (x is row vector),
    % a' = Da (a is column vector) etc.

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.

    % Check inputs
    '''
    assert isinstance(nmax, float) or isinstance(nmax, int), 'Nmax must be a numeric scalar'
    assert isinstance(R, np.ndarray) and R.shape == (3,3), 'R must be a 3x3 rotation matrix'
    C = np.array([[1/np.sqrt(2), 0, -1/np.sqrt(2)],
        [1j/np.sqrt(2), 0, 1j/np.sqrt(2)],
        [0, 1, 0]])
    invC = np.array([[1/np.sqrt(2), -1j/np.sqrt(2), 0],
                    [0, 0, 1],
                    [-1/np.sqrt(2), -1j/np.sqrt(2), 0]])
    D = np.matmul(invC, np.matmul(R, C))
    D1 = D.T
    DD = D1.copy()
    X = [D1]
    for i in range(2, nmax+1):
        DDD = 0j + np.zeros((2*i+1, 2*i+1))
        m0 = np.matmul(np.ones((2*i-1,1)), np.arange(-i, i+1).reshape((1, -1)))
        m1 = np.arange(-i+1, i).reshape((-1, 1))*np.ones((1, 2*i+1))
        a = np.sqrt((i + m0)*(i - m0)/((i + m1)*(i - m1)))
        b = np.sqrt((i + m0)*(i + m0 - 1)/( 2*(i + m1)*(i - m1)))
        DDD[1:-1, 1:-1] = D1[1, 1] * a[:, 1:-1]*DD
        DDD[1:-1, 2:] = DDD[1:-1, 2:] + D1[1,2] * b[:,2:]*DD
        DDD[1:-1, :-2] = DDD[1:-1,:-2] + D1[1,0]*np.fliplr(b[:, 2:])*DD
        m0 = np.arange(-i, i+1).reshape((1, -1))
        c = np.sqrt((i + m0)*(i - m0)/(i*(2*i - 1)))
        d = np.sqrt((i + m0)*(i + m0 - 1)/(2*i*(2*i - 1)))
        DDD[0, 1:-1] = D1[0, 1]*c[0, 1:-1]*DD[0,:]
        DDD[0, 2:] = DDD[0,2:] + D1[0,2]*d[:, 2:]*DD[0,:]
        DDD[0, :-2] = DDD[0, :-2] + D1[0, 0] * np.fliplr(d[:, 2:])*DD[0, :]  
        DDD[-1, 1:-1] = D1[2, 1] * c[:, 1:-1] * DD[-1, :]
        DDD[-1, 2:] = DDD[-1, 2:] + D1[2,2] * d[:, 2:]*DD[-1, :]
        DDD[-1,:-2] = DDD[-1,:-2] + D1[2,0] * np.fliplr(d[:, 2:]) * DD[-1, :]
        DD = DDD.copy()
        X.append(DDD)        
    return block_diag(*X)
