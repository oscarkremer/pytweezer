import numpy as np

def combined_index(in1, in2=None):
    '''
        %COMBINED_INDEX translates between (n,m) and combined index
    % Mode indices and combined index are related by: ci = n * (n+1) + m.
    %
    % [n,m] = COMBINED_INDEX(ci) calculates (n,m) from the combined index.
    %
    % ci = COMBINED_INDEX(n,m) calculates the combined index from mode indices.
    %
    % length = COMBINED_INDEX(Nmax, Nmax) calculates length of the beam vectors.
    %
    % Nmax = COMBINED_INDEX(length) calculates Nmax from length of beam vectors.

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.

    % Sanity check
    '''
    if isinstance(in2, np.ndarray): 
        return in1*(in1 + 1) + in2
    else:
        if not in2:
            out1 = np.floor(np.sqrt(in1))
            out2 = in1 - np.power(out1, 2) - out1
            return out1, out2
        else:
            return in1*(in1 + 1) + in2    
        