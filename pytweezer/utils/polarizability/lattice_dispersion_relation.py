import numpy as np
from .clausius_mossoti import clausius_mossoti

def lattice_dispersion_relation(spacing, index, kvec=np.array([]), E=np.array([]), k=2*np.pi):
    '''  % Lattice dispersion relation polarizablity
    %
    % Polarizability calculation based on
    %
    %   Draine & Goodman, Beyond Clausius-Mossoti: wave propagation
    %   on a polarizable point lattice and the discrete dipole approximation,
    %   The Astrophysical Journal, 405:685-697, 1993 March 10
    %
    % Usage
    %   alpha = LDR(spacing, index, ...)
    %   Calculates a Nx1 element vector containing the isotropic
    %   polarisabilities for N dipoles.
    %
    %   alpha = LDR(spacing, index, kvec, E0, ...)
    %   As above but specifies the polarisability information for use
    %   with plane wave illumination.
    %
    % Parameters
    %   - spacing (numeric scalar) -- lattice spacing parameter
    %   - index (Nx1 numeric) -- Relative refractive indices for N dipoles.
    %   - kvec (1x3 numeric) -- Wave vector [kx, ky, kz]
    %   - E0 (1x3 numeric) -- E-field polarisation [Ex, Ey, Ez]
    %
    % Optional named arguments
    %   - k0 (numeric) -- Wavenumber to scale spacing by.  Default: ``2*pi``.

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
            b = np.array([-1.8915316, 0.1648469, -1.7700004])
            if kvec.size and E.size:
                if kvec.shape != E.shape:
                    raise ValueError('Mismatch between dimensions of E and kvec variables')
                else:
                    a_hat = np.power(kvec/np.linalg.norm(kvec), 2)
                    e_hat = np.power(E/np.linalg.norm(E), 2)
                    S = np.matmul(a_hat.T, e_hat)
            else:
                S = 0.2
            alpha_cm = clausius_mossoti(spacing, index)
            alpha = alpha_cm/(1 + (alpha_cm/dcube)*((b[0]+msqr*b[1]+msqr*b[2]*S)*(k*spacing)**2-2/3*1j*k**3*dcube))
            return alpha
        else:
            raise ValueError('Index variable must be an 1 or 2 dimensional array')
    else:
        raise TypeError('Index variable must be numpy array')
