import numpy as np

def clausius_mossoti(spacing: float, index: np.ndarray):
    '''    % Clausius-Mossoti Polarizability 
    %
    % Usage
    %   alpha = CM(spacing, index)
    %   Calculates a Nx1 element vector containing the isotropic
    %   polarisabilities for N dipoles.
    %
    % Parameters
    %   - spacing (numeric scalar) -- lattice spacing parameter
    %   - index (Nx1 numeric) -- Relative refractive indices for N dipoles.

    % Based on the script by Vincent Loke.
    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.
    '''
    if isinstance(index, np.ndarray):
        index_shape = index.shape
        if len(index_shape) in [1, 2]:
            index = index.reshape((max(index_shape), 1))
            msqr = np.power(index ,2)
            dcube = spacing**3
            alpha = (3*dcube/(4*np.pi))*(msqr - 1)/(msqr + 2)
            return alpha        
        else:
            raise ValueError('Index variable must be an 1 or 2 dimensional array')
    else:
        raise TypeError('Index variable must be numpy array')
