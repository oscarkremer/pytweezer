def paraxial_beam_waist(paraxial_order, tolerance=0.00001):
    '''% PARAXIAL_BEAM_WAIST computes the re-normlaised beam waist for high-order
    % gaussian beams using a recursion to find a particular root of the function
    % x.^\lambda exp(-x^2).
    %
    % w = paraxial_beam_waist(paraxial_order) computes the beam waist parameter
    % for input into a paraxial code.

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.
    '''
    w = 1.0
    
    if paraxial_order:
        inv_L = 1/abs(paraxial_order)
        zz = np.exp(-(abs(paraxial_order)+2)*invL)
        w = -(1 + 2*np.sqrt(inv_L)+inv_L);
        
        w0 = -w
        
        while (abs(w-w0) > tolerance):
            w0 = w
            expw = np.exp(w)  
            w=w0-(w0*expw+zz)/(expw+w0*expw); %Newton's rule... Usually this kind of method would find the real root i.e. W_0(z)... This finds W_{-1}(z) local to the beam waist of an LG beam.
        
        w = sqrt(-abs(paraxial_order )/2.*w); %Beam waist in normalized units
    return w
