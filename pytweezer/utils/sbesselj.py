import numpy as np
from scipy.special import jv

def sbesselj(n, kr, verbose=False):
    '''    % SBESSELJ spherical bessel function jn(kr)
    % jn(kr) = sqrt(pi/2kr) Jn+0.5(kr)
    %
    % jn = SBESSEL(n,z) calculates the spherical bessel function.
    %
    % [jn,dzjn] = sbessel(n,z) additionally, calculates the derivative
    % of the appropriate Ricatti-Bessel function divided by z.
    %
    % See also besselj.

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.
    '''

    if isinstance(n, np.ndarray):
        n = np.concatenate([n, n-1])
    elif isinstance(n, float) or isinstance(n, int):
        n = np.array([n, n-1])
    else:
        raise TypeError('Parameter \'n\' must numpy.ndarray, integer or float!')

    n, kr = np.meshgrid(n, kr)
    jn = jv(n+1/2, kr)
    small_args = np.argwhere(np.abs(kr) < 1e-15)
    not_small_args = np.argwhere(np.abs(kr) > 1e-15)
    if (isinstance(kr, float) or isinstance(kr, int)) and abs(kr) < 1e-15:
        if verbose:
            print('Computing spherical bessel function for single parameter and small value!')
        jn = kr**n / np.prod(np.arange(1, 2*n+2, 2))
    elif (isinstance(kr, float) or isinstance(kr, int)) and abs(kr) > 1e-15:
        if verbose:
            print('Computing spherical bessel function for single parameter and non-small value!')
        jn = np.sqrt(np.pi/(2*kr))* jn
    elif (isinstance(n, float) or isinstance(n, int)) == 1:
        jn[not_small_args] = np.sqrt(np.pi/(2*kr(not_small_args)))* jn[not_small_args]
        jn[small_args] = kr[small_args]^n / np.prod(np.arange(1, 2*n+2, 2))
    else:
        if not_small_args.size:
            not_small_args = tuple(not_small_args.T)    
            jn[not_small_args] = np.sqrt(np.pi/(2*kr[not_small_args]))*jn[not_small_args]
        if small_args.size:
            small_args = tuple(small_args.T)

            jn[small_args] = kr[small_args]**n[small_args].astype(float)/np.prod(np.arange(1, 2*n[small_args].all()+2, 2))
    cols = int(jn.shape[1]/2)
    djn = jn[:, cols:] - n[:, :cols]/kr[:, :cols]*jn[:, :cols]

    jn = jn[:, :cols]
    return jn, djn
