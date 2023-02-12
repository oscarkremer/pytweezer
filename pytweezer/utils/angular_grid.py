import numpy as np


def angular_grid(ntheta, nphi, style='column'):
    '''
    %ANGULARGRID makes a angular grid of points over a sphere
    %
    % [theta,phi] = ANGULARGRID(N) generates two column N^2-by-1 matrices
    % with theta (polar) and phi (azimuthal) angle pairs for N discrete
    % evenly spaced polar and azimuthal angles.
    %
    % [theta,phi] = ANGULARGRID(ntheta, nphi) specifies the number of
    % evenly spaced points to use in the theta and phi direction.
    %
    % [theta,phi] = ANGULARGRID(..., behaviour) uses behaviour to control
    % the output type:
    %
    % behaviour | output
    % -----------------------------------------------
    %     0     | column vectors of all points
    %     1     | vectors of all theta and phi values 
    %     2     | ntheta x nphi matrix of all points
    %
    % Note that the output data values are the same for
    % behaviours 0 and 2; they're just arranged differently.
    % To convert from one format to another:
    % 2 -> 0: theta = theta(:); phi = phi(:);
    % 0 -> 2: theta = reshape(theta,ntheta,nphi);
    %         phi = reshape(phi,ntheta,nphi);

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.
    '''
    # theta goes from 0 to pi - we avoid the endpoints
    # since they are mathematically troublesome
    if style not in ('column', 'points', 'matrix'):
        raise ValueError('Style for grid creation not allowed!!! Allowed values are: (`column`, `points` and `matrix`)')
    else:
        theta = np.pi*(np.arange(1, ntheta+1)-0.5)/ntheta
        theta = theta.reshape((ntheta, 1)).T
        phi = 2*np.pi*(np.arange(1, nphi+1)-1)/nphi
        phi = phi.reshape((nphi, 1)).T
        unitary_matrix = np.ones((nphi, ntheta))
        theta = (unitary_matrix*theta).T
        phi = (unitary_matrix.T*phi)
        if style == 'matrix':
            return theta, phi
        elif style == 'points':
            return theta.T.reshape((theta.shape[0]*theta.shape[1], 1)), phi.T.reshape((phi.shape[0]*phi.shape[1], 1))     
        else:
            return theta.T.flatten(), phi.T.flatten()
