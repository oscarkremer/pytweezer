import numpy as np

def legendre_row(n, theta):
    '''LEGENDREROW gives the spherical coordinate recursion in m
    %
    % pnm = LEGENDREROW(n, theta) gives the spherical recursion for a given
    % n, theta.
    %
    % This provides approximately no benefit over the MATLAB implimentation. It
    % *may* provide a benefit in Octave. Inspiration from
    % [Holmes and Featherstone, 2002] and [Jekeli et al., 2007].

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.
    '''
    if n==0:
        pnm=1/np.sqrt(2*pi)/np.sqrt(2)
        return pnm


    ct = np.cos(theta)
    st = np.sin(theta)

    Wnn = np.sqrt((2*n+1)/(4*np.pi)*np.prod(1-1/2/np.arange(1,n+1)))*np.ones(theta.shape)
    Wnnm1=np.sqrt(2*n)*ct*Wnn 
    lnm = np.arange(0, n+1).shape[0]
    pnm = np.zeros((lnm, theta.shape[0]))
    pnm[-1,:] = Wnn
    pnm[-2,:] = Wnnm1
    if lnm == 2:
        pass
    else:
        jj = lnm-2
        for i in range(n-2,-1,-1):
            a = np.sqrt(4*(i+1)**2/(n-i)/(n+i+1))
            b = np.sqrt((n-i-1)*(n+i+2)/(n-i)/(n+i+1))


            pnm[jj-1,:] = a*ct*pnm[jj,:]-b*st**2*pnm[jj+1,:]
            jj = jj-1
    ST, M = np.meshgrid(st, np.arange(0, n+1))
    
    pnm = pnm*ST**M
    return pnm
  