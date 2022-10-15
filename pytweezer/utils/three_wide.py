import numpy as np


def three_wide(a):
    '''    % THREEWIDE creates colum vector with input repeated in 3 columns
    %   the function can take a column of row vector input, the output
    %   will be a matrix with three columns.
    %
    % You might find this useful for multiplying a vector of scalars
    % with a column vector of 3-vectors.

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.
    '''
    if isinstance(a, float) or isinstance(a, int):
        return np.array([[a, a, a]])
    elif isinstance(a, np.ndarray):
        if len(a.shape) == 1:
            return np.array([a, a a]).T 
        else:
            raise ValueError('Numpy array must have one dimensional ')
    else:
        raise ValueError('Input must be one dimensional numpy array, integer or float.')