import numpy as np
from .match_size import match_size
from .spherical_harmonics import spherical_harmonics

def vsh(n, m, theta, phi):
    '''
        % VSH calculate vector spherical harmonics
        %
        % [B,C,P] = VSH(n,m,theta,phi) calculates vector spherical harmonics
        % for the locations theta, phi.  Vector m allowed.  Scalar n for the moment.
        %
        % [B,C,P] = VSH(n,theta,phi) outputs for all possible m.
        %
        % If scalar m: B,C,P are arrays of size length(theta,phi) x 3
        % If vector m: B,C,P are arrays of size length((theta,phi),m) x 3
        % theta and phi can be vectors (of equal length) or scalar.
        %
        % The three components of each vector are [r,theta,phi]
        %
        % "Out of range" n and m result in return of [0 0 0]

        % This file is part of the optical tweezers toolbox.
        % See LICENSE.md for information about using/distributing this file.
    '''

    if not isinstance(n, int) and not isinstance(n, float):
        raise TypeError('Input parameter \'n\' must be scalar.')
    theta, phi = match_size(theta, phi)
    Y, Ytheta, Yphi = spherical_harmonics(n, m, theta, phi)
    Z = np.zeros(Y.shape)
    B = np.concatenate([Z, Ytheta, Yphi])
    C = np.concatenate([Z, Yphi, -Ytheta])
    P = np.concatenate([Y, Z, Z])
    return B, C, P