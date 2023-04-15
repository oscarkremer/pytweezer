import numpy as np
from .meshgrid import meshgrid
from ._nb_legendre_row_ import _nb_legendre_row_

def legendre_row(n, theta):
    '''LEGENDREROW gives the spherical coordinate recursion in m
    %
    % pnm = LEGENDREROW(n, theta) gives the spherical recursion for a given
    % n, theta.
    %
    % This provides approximately no benefit over the MATLAB implimentation. It
    % *may* provide a benefit in Octave. Inspiration from
    % [Holmes and Featherstone, 2002] and [Jekeli et al., 2007].

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.
    '''
    st, pnm = _nb_legendre_row_(n, theta)
    ST, M = meshgrid(st, np.arange(0, n+1))
    pnm = pnm*np.power(ST, M)
    return pnm
  