import numpy as np

def laguerre(p, l, x):
    '''% LAGUERRE associated Laguerre function
    %
    % L = LAGUERRE(p,l,X) evaluate associated Laguerre function.
    % p and l must be integer scalars greater than zero

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.
    '''
    l_plt = np.zeros((x.size, max(p+1,2)))
    l_plt[:, 0] = 1
    l_plt[:, 1] =  l - x + 1
    for i in range(2, p+1):
        l_plt[:,i] = 1/i*((2*i+l-1-x)*Lplt[:, i-1]-(i+l-1)*l_plt[:, i-2])
    l_pl = l_plt[:, p].reshape(x.shape)
    return l_pl