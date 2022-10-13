from scipy.special import hankel1, hankel2


def spherical_hankel(kr, n, htype):
    '''    % SBESSELH spherical hankel function hn(kr) of the first or 
    % second kind hn(kr) = sqrt(pi/2kr) Hn+0.5(kr)
    %
    % hn = SBESSELH(n,htype,z) computes the spherical hankel function
    % of degree, n, of kind, htype, and argument, z. 
    %
    % [hn,dzhn] = SBESSELH(n,htype,z) additionally, calculates the 
    % derivative of the appropriate Ricatti-Bessel function divided 
    % by z.
    %
    % See also besselj and bessely.

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.'''

    #ott.warning('internal');
    if htype not in ('1', '2'):
        raise ValueError('Type of Hankel function must be `1` or `2`')
    kr=kr.T
    n=n.T

    n, kr=np.meshgrid(n, kr)
    if htype='1':
        hn = hankel1(n+1/2, kr)
    elif htype='2':
        hn = hankel2(n+1/2, kr)    
    hn = np.sqrt(np.pi/(2*kr)) * (hn)


    dhn = hn[:, int((hn.shape[1]-1)/2)+1:]-n[:,:int(n.shape[1]/2)]/kr[1:,:int(n.shape[1]/2)]*hn[:,1:int(n.shape[1]/2)]
    hn = hn[:,:int((hn.shape[1]-1)/2)]
    
    #ott.warning('external');
    return hn,dhn