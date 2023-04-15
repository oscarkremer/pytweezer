import numpy as np
import warnings
warnings.filterwarnings("ignore", category=PendingDeprecationWarning)
from numpy import matlib as matlib
from .match_size import match_size
from .sbesselh import sbesselh
from .sbesselj import sbesselj
from .three_wide import three_wide
from .vsh import vsh
from numba import jit

#@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
def vswf(n, m, kr, theta, phi, htype=None):
    '''
    % VSWF vector spherical wavefunctions: M_k, N_k.
    %
    % [M1,N1,M2,N2,M3,N3] = VSWF(n,m,kr,theta,phi) calculates the
    % outgoing M1,N1, incomming M2,N2 and regular M3,N3 VSWF.
    % kr, theta, phi are vectors of equal length, or scalar.
    %
    % [M,N] = VSWF(n,m,kr,theta,phi,type) calculates only the
    % requested VSWF, where type is
    %     1 -> outgoing solution - h(1)
    %     2 -> incoming solution - h(2)
    %     3 -> regular solution - j (ie RgM, RgN)
    %
    % VSWF(n, kr, theta, phi) if m is omitted, will calculate for all m.
    %
    % M,N are arrays of size length(vector_input,m) x 3
    %
    % The three components of each vector are [r,theta,phi].
    %
    % "Out of range" n and m result in return of [0 0 0]

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.
    '''
    #[M,N,M2,N2,M3,N3] 

    if not (isinstance(n, float) or isinstance(n, int)):
        raise TypeError('Variable \'n\' must integer of float.')
    htype = 0 if not htype else htype
    if isinstance(htype, str):        
        if htype == 'incoming':
            htype = 2
        elif htype == 'outgoing':
            htype = 1
        elif htype == 'regular':
            htype = 3
        else:
            raise ValueError('Unknown htype string')


    kr, theta, phi = match_size(kr, theta, phi)

    B, C, P = vsh(n, m, theta, phi)
    if n > 0:
        Nn = np.sqrt(1/(n*(n+1)))
    else:
        Nn = 0
    if htype == 1:
        if not (isinstance(m, int) or isinstance(m, float) or isinstance(m, np.int64)):
            kr3 = matlib.repmat(kr,[1, m.shape[0]*3])
            hn = matlib.repmat(sbesselh(n, kr, '1')[0],[1, m.shape[0]*3])
            hn1 = matlib.repmat(sbesselh(n-1, kr, '1')[0],[1, m.shape[0]*3])
        else:
            kr3 = three_wide(kr)
            hn = three_wide(sbesselh(n, kr, '1')[0])
            hn1 = three_wide(sbesselh(n-1, kr, '1')[0])
        M = Nn*hn*C
        N = Nn*(n*(n+1)/kr3*hn*P+(hn1-n/kr3*hn)*B)
        M2 = 0
        N2 = 0
        M3 = 0
        N3 = 0
    elif htype == 2:
        if not (isinstance(m, int) or isinstance(m, float) or isinstance(m, np.int64)):
            kr3 = matlib.repmat(kr,[1, m.shape[0]*3])
            hn = matlib.repmat(sbesselh(n,kr, '2')[0], [1, m.shape[0]*3])
            hn1 = matlib.repmat(sbesselh(n-1,kr, '2')[0], [1, m.shape[0]*3])
        else:
            kr3 = three_wide(kr)
            hn = three_wide(sbesselh(n, kr, '2')[0])
            hn1 = three_wide(sbesselh(n-1,kr, '2')[0])
            
        M = Nn * hn * C
        N = Nn * ( n*(n+1)/kr3*hn * P + (hn1 - n/kr3*hn)*B)
        M2 = 0
        N2 = 0
        M3 = 0
        N3 = 0
    elif htype == 3:
        if not (isinstance(m, int) or isinstance(m, float) or isinstance(m, np.int64)):
            kr3 = matlib.repmat(kr,[1, m.shape[0]*3])
            jn = matlib.repmat(sbesselj(n, kr)[0],[1, m.shape[0]*3])
            jn1 = matlib.repmat(sbesselj(n-1, kr)[0],[1, m.shape[0]*3])
        else:
            kr3 = three_wide(kr)
            jn = three_wide(sbesselj(n, kr)[0])
            jn1 = three_wide(sbesselj(n-1, kr)[0])

        M = Nn * jn * C
        N = Nn * (n*(n+1)/kr3*jn*P + ( jn1 - n/kr3*jn)* B)
        M2 = 0
        N2 = 0
        M3 = 0
        N3 = 0
        if n != 1: 
            N[kr3==0]=0
        else:
            N[kr3==0] = 2/3*Nn*( P.reshape((1, -1))[kr3==0] + B.reshape((1,-1))[kr3==0])
    else:
        if not (isinstance(m, int) or isinstance(m, float) or isinstance(m, np.int64)):
            kr3 = matlib.repmat(kr,[1, m.shape[0]*3])
            jn = matlib.repmat(sbesselj(n, kr)[0],[1, m.shape[0]*3])
            jn1 = matlib.repmat(sbesselj(n-1, kr)[0],[1, m.shape[0]*3])
            hn1 = matlib.repmat(sbesselh(n, kr, '1')[0],[1, m.shape[0]*3])
            hn11 = matlib.repmat(sbesselh(n-1, kr, '1')[0], [1, m.shape[0]*3])            
            hn2 = matlib.repmat(sbesselh(n, kr, '2')[0],[1, m.shape[0]*3])
            hn21 = matlib.repmat(sbesselh(n-1, kr, '2')[0],[1, m.shape[0]*3])
        else:
            kr3 = three_wide(kr)
            hn2 = three_wide(sbesselh(n, kr, '2')[0])
            hn21 = three_wide(sbesselh(n-1, kr, '2')[0])  
            hn1 = three_wide(sbesselh(n, kr, '1')[0])
            hn11 = three_wide(sbesselh(n-1, kr, '1')[0]) 
            jn = three_wide(sbesselj(n, kr)[0])
            jn1 = three_wide(sbesselj(n-1, kr)[0])
        
        M = Nn * hn1 * C
        N = Nn * ( n*(n+1)/kr3* hn1 * P + ( hn11 - n/kr3 * hn1 ) * B )
        M2 = Nn * hn2 * C
        N2 = Nn * ( n*(n+1)/kr3 * hn2 * P + ( hn21 - n/kr3 * hn2 ) * B )
        M3 = Nn * jn* C
        N3 = Nn * ( n*(n+1)/kr3 * jn * P + ( jn1 - n/kr3 * jn ) * B )
        if n!=1:
            N3[kr3==0] = 0
        else:
            N3[kr3==0] = 2/3*Nn*( P[kr3==0] + B[kr3==0])

    return M, N, M2, N2, M3, N3