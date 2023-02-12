def n_max_2ka(n_max):
    '''    % NMAX2KA finds size parameter ka corresponding to Nmax
    %
    % ka = NMAX2KA(Nmax) finds size parameter for maximum order, Nmax,
    % which spherical expansions are truncated.
    %
    % Truncation order is given by Nmax = ka + 3 (ka)^(1/3)

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.
    '''
    for i in range(1, length(n_max)):
        kas = roots([1, (-3*n_max[i]), (27+3*n_max[i]**2), (-n_max[i]**3)]);
        ka[i] = kas[2]
    return ka 
