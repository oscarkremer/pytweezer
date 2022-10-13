import numpy as np
from .clausius_mossoti import clausius_mossoti

def filter_coupled_dipole(spacing: float, index: np.ndarray,
                        k=2*np.pi):
    '''
    % Filtered coupled dipole polarizability
    %
    % Usage
    %   alpha = FCD(spacing, index)
    %   Calculates a Nx1 element vector containing the isotropic
    %   polarizabilities for N dipoles.
    %
    % Parameters
    %   - spacing (numeric scalar) -- lattice spacing parameter
    %   - index (Nx1 numeric) -- Relative refractive indices for N dipoles.
    %
    % Optional named arguments
    %   - k0 (numeric) -- Wavenumber to scale spacing by.  Default: ``2*pi``.

    % Based on the script by Vincent Loke.
    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.

    p = inputParser;
    p.addParameter('k0', 2*pi);
    p.parse(varargin{:});

    k = p.Results.k0;
    msqr = index(:).^2;

    alpha_CM = 3*spacing^3/(4*pi)*(msqr - 1)./(msqr + 2); % Clausius-Mossotti
        alpha_FCD = alpha_CM./(1 + (alpha_CM/spacing^3).*(4/3*(k*spacing)^2 + 2/3*(1i + log((pi-k*spacing)/(pi+k*spacing))/pi)*k^3*spacing^3));
    return alpha_FCD
    '''
    if isinstance(index, np.ndarray):
        index_shape = index.shape
        if len(index_shape) in [1, 2]:
            alpha_cm = clausius_mossoti(spacing, index)
            alpha = alpha_cm/(1+(alpha_cm/spacing**3)*(4/3*(k*spacing)**2 + 2/3*(1j + np.log((0j+np.pi-k*spacing)/(0j+np.pi+k*spacing))/np.pi)*((k*spacing)**3)))
            return alpha
        else:
            raise ValueError('Index variable must be an 1 or 2 dimensional array')
    else:
        raise TypeError('Index variable must be numpy array')