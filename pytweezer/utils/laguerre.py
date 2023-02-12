def laguerre(p, l, X):
    '''% LAGUERRE associated Laguerre function
    %
    % L = LAGUERRE(p,l,X) evaluate associated Laguerre function.
    % p and l must be integer scalars greater than zero

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.
    '''
    Lplt = np.zeros((x.size, max(p+1,2)))
    Lplt[:, 0] = 1
    Lplt[:, 1] =  l - x + 1
    for i in range(2, p+1):
        Lplt[:,i] = 1/i*((2*i+l-1-x)*Lplt[:, i-1]-(i+l-1)*Lplt[:, i-2])
    Lpl = Lplt[:, p].reshape(X.shape)
    return Lpl