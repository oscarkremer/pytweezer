import numpy as np
from .rtp2xyz import rtp2xyz
from .vswf import vswf
from numba import jit

#@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
def vswf_cart(n, m, kr, theta, phi, htype=None):
    '''
    % VSWFCART vector spherical harmonics spherical coordinate input,
    % cartesian output.
    %
    % [M1,N1,M2,N2,M3,N3] = VSWFCART(n,m,kr,theta,phi) calculates the
    % outgoing M1,N1, incomming M2,N2 and regular M3,N3 VSWF.
    % kr, theta, phi are vectors of equal length, or scalar.
    %
    % [M,N] = VSWFCART(n,m,kr,theta,phi,type) calculates only the
    % requested VSWF, where type is
    %     1 -> outgoing solution - h(1)
    %     2 -> incoming solution - h(2)
    %     3 -> regular solution - j (ie RgM, RgN)
    %
    % Scalar n,m for the moment.
    % M,N are arrays of size length(vector_input) x 3
    %
    % The three components of each input vector are [kr,theta,phi]
    % The three components of each output vector are [x,y,z]
    %
    % "Out of range" n and m result in return of [0 0 0]
    %
    % At the coordinate origin (kr == 0) we use only theta/phi.

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.
    '''

    if not htype:
        htype = 0
    
    M, N, _, _, _, _ = vswf(n, m, kr, theta, phi, htype)
    x, y, z = rtp2xyz(kr,theta,phi)
    theta_hat_x = np.cos(theta)*np.cos(phi)
    theta_hat_y = np.cos(theta)*np.sin(phi)
    theta_hat_z = -np.sin(theta)
    phi_hat_x = -np.sin(phi)
    phi_hat_y = np.cos(phi)
    phi_hat_z = 0

    kr_safe = kr
    if not isinstance(kr, np.ndarray):
        if abs(kr) < 1e-15:
            kr_safe = 1.0
    else:
        kr_safe[kr==0] = 1.0
    r_hat_x = x/kr_safe
    r_hat_y = y/kr_safe
    r_hat_z = z/kr_safe
    if not isinstance(kr, np.ndarray):
        if abs(kr) < 1e-15:
            r_hat_x = np.sin(theta)*np.cos(phi)
            r_hat_y = np.sin(theta)*np.sin(phi)
            r_hat_z = np.cos(theta)
    else:
        r_hat_x[kr==0] = np.sin(theta[kr == 0])*np.cos(phi[kr==0])
        r_hat_y[kr==0] = np.sin(theta[kr == 0])*np.sin(phi[kr==0])
        r_hat_z[kr==0] = np.cos(theta[kr == 0])
    outputs = []
    for matrix in [M, N]:
    
        npts = int(matrix.shape[1]/3)
        Mr = matrix[:, :int(npts)]
        Mtheta = matrix[:,npts:2*npts]
        Mphi = matrix[:, 2*npts:3*npts]
        Mx = Mr* r_hat_x + Mtheta* theta_hat_x + Mphi * phi_hat_x
        My = Mr* r_hat_y + Mtheta* theta_hat_y + Mphi * phi_hat_y
        Mz = Mr* r_hat_z + Mtheta* theta_hat_z + Mphi * phi_hat_z
        outputs.append(np.array([Mx, My, Mz]))
    return outputs[0], outputs[1]