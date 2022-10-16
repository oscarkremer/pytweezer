import numpy as np
from scipy.special import hankel1, hankel2


def sbesselh(n, kr, htype):
    '''    % SBESSELH spherical hankel function hn(kr) of the first or 
    % second kind hn(kr) = sqrt(pi/2kr) Hn+0.5(kr)
    %
    % hn = SBESSELH(n,htype,z) computes the spherical hankel function
    % of degree, n, of kind, htype, and argument, z. 
    %
    % [hn,dzhn] = SBESSELH(n,htype,z) additionally, calculates the 
    % derivative of the appropriate Ricatti-Bessel function divided 
    % by z.
    %
    % See also besselj and bessely.

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.'''

    #ott.warning('internal');
    if htype not in ('1', '2'):
        raise ValueError('Type of Hankel function must be `1` or `2`')
    if isinstance(n, np.ndarray):
        n = np.concatenate([n, n-1])
    elif isinstance(n, float) or isinstance(n, int):
        n = np.array([n, n-1])
    else:
        raise TypeError('Parameter \'n\' must numpy.ndarray, integer or float!')

    n, kr = np.meshgrid(n, kr)
    if htype == '1':
        hn = hankel1(n+1/2, kr)
    elif htype == '2':
        hn = hankel2(n+1/2, kr)    
    hn = np.sqrt(np.pi/(2*kr)) * hn
    cols = int(hn.shape[1]/2)
    dhn = hn[:, cols:]-n[:, :cols]/kr[:, :cols]*hn[:, :cols]
    hn = hn[:, :cols]
    return hn, dhn