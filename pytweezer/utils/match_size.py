import warnings
warnings.filterwarnings("ignore", category=PendingDeprecationWarning)
import numpy as np
from numpy import matlib as matlib


def match_size(*args):
    '''    % Checks that all vector inputs have the same number of rows.
    %
    % Usage
    %   [A,B,...] = matchsize(A,B,...) checks inputs have same number of rows,
    %   and expands single-row inputs by repetition to match the input row number.
    %
    % Parameters
    %   - A,B,... (numeric)  -- Numeric arrays whose number of rows
    %     are to be matched.
    %
    % Example
    %   The following example shows has two inputs, a scalar and a row vector.
    %   The scalar is duplicated to match the length of the row vector::
    %
    %     A = 5;
    %     B = [1; 2; 3];
    %
    %     [A,B] = matchsize(A, B);
    %     disp(A)  % -> [5; 5; 5]
    %     disp(B)  % -> [1; 2; 3]
    '''
    #% This file is part of the optical tweezers toolbox.
    #% See LICENSE.md for information about using/distributing this file.

    #% Loop over inputs to determine maximum size
    nmax = 0
    for arg in args: 
        if isinstance(arg, np.int64) or isinstance(arg, int) or isinstance(arg, float):
            nmax = max(1, nmax)
        elif arg.shape[0] > nmax:
            nmax = max(nmax, arg.shape[0])
  
    for arg in args:
        if isinstance(arg, int) or isinstance(arg, float) or isinstance(arg, np.int64):
            nrows = 1
        else:
            nrows = arg.shape[0]
        if nrows == 1:
           yield matlib.repmat(arg, nmax, 1)
        elif nrows != nmax:
            raise ValueError('Number of rows in inputs must be one or equal.')
        else:
            yield arg