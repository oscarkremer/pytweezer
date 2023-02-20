import numpy as np
from .legendre_row import legendre_row
from .match_size import match_size

def spherical_harmonics(n, m, theta, phi=np.array([])):
    '''
    % SPHARM scalar spherical harmonics and angular partial derivatives.
    %
    % Y = SPHARM(n,m,theta,phi) calculates scalar spherical harmonics.
    %
    % [Y,Ytheta,Yphi] = SPHARM(n,m,theta,phi) additionally, calculates
    % the angular partial derivatives dY/dtheta and 1/sin(theta)*dY/dphi.
    %
    % SPHARM(n,theta,phi) as above but for all m.
    %
    % Scalar n for the moment.
    % 
    % If scalar m is used Y is a vector of length(theta,phi) and is
    % completely compatible with previous versions of the toolbox. If vector m
    % is present the output will be a matrix with rows of length(theta,phi) for
    % m columns.
    %
    % "Out of range" n and m result in return of Y = 0

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.
    '''
    if not (isinstance(n, float) or isinstance(n, int)):
        raise TypeError('Input parameter \'n\' must be scalar.')
    if not phi.size:    
        phi = theta
        theta = m
    mi = m
    m = np.arange(-n,n+1)
    index_bigger = np.argwhere(abs(m) <=n)
    if index_bigger.size:
       m = m[tuple(index_bigger.T)]
    theta, phi = match_size(theta, phi)
    input_length = theta.shape[0]
    pnm = legendre_row(n, theta)

    pnm = pnm[abs(m), :]
    
    phiM, mv = np.meshgrid(phi, m)
    index_m_positive = np.argwhere(m>=0)
    index_m_negative = np.argwhere(m<0)

    pnm = np.concatenate([1/((-1)**(-mv[tuple(index_m_negative),:]))*pnm[tuple(index_m_negative),:], pnm[tuple(index_m_positive),:]])
    print(pnm.shape)

    pnm = pnm.reshape((pnm.shape[0], -1))
    expphi = np.exp(1j*mv*phiM)
    Y = pnm*expphi
    expplus = np.exp(1j*phiM)
    expminus = np.exp(-1j*phiM)

    ymplus = np.concatenate([Y[1:,:], np.zeros((1,theta.shape[0]))])
    ymminus = np.concatenate([np.zeros((1,theta.shape[0])), Y[:-1,:]])
    Ytheta = np.sqrt((n-mv+1)*(n+mv))/2*expplus*ymminus - np.sqrt((n-mv)*(n+mv+1))/2*expminus*ymplus

    Y2 = _spherical_harmonics(n+1,theta,phi)
    ymplus = Y2[2:, :]
    ymminus = Y2[:-2,:]
    Yphi = 1j/2* np.sqrt((2*n+1)/(2*n+3))*(np.sqrt((n+mv+1)*(n+mv+2))*expminus*ymplus+np.sqrt((n-mv+1)*(n-mv+2))*expplus*ymminus)
    Y = Y[n+mi,:]
    Yphi = Yphi[n+mi,:]
    Ytheta = Ytheta[n+mi,:]
    return Y, Ytheta, Yphi



def _spherical_harmonics(n, m, theta, phi=None):
    '''
    % SPHARM scalar spherical harmonics and angular partial derivatives.
    %
    % Y = SPHARM(n,m,theta,phi) calculates scalar spherical harmonics.
    %
    % [Y,Ytheta,Yphi] = SPHARM(n,m,theta,phi) additionally, calculates
    % the angular partial derivatives dY/dtheta and 1/sin(theta)*dY/dphi.
    %
    % SPHARM(n,theta,phi) as above but for all m.
    %
    % Scalar n for the moment.
    % 
    % If scalar m is used Y is a vector of length(theta,phi) and is
    % completely compatible with previous versions of the toolbox. If vector m
    % is present the output will be a matrix with rows of length(theta,phi) for
    % m columns.
    %
    % "Out of range" n and m result in return of Y = 0

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.
    '''
    if not (isinstance(n, float) or isinstance(n, int)):
        raise TypeError('Input parameter \'n\' must be scalar.')
    if not phi:    
        phi = theta
        theta = m
    mi = m
    m = np.arange(-n,n+1)
    index_bigger = np.argwhere(abs(m) <=n)
    if index_bigger.size:
       m = m[tuple(index_bigger.T)]
    theta, phi = match_size(theta, phi)
    input_length = theta.shape[0]
    pnm = legendre_row(n, theta)
    pnm = pnm[abs(m), :] 
    phiM, mv = np.meshgrid(phi, m)
    index_m_positive = np.argwhere(m>=0)
    index_m_negative = np.argwhere(m<0)
    pnm = np.concatenate([1/((-1)**(-mv[tuple(index_m_negative),:]))*pnm[tuple(index_m_negative),:], pnm[tuple(index_m_positive),:]])
    
    pnm = pnm.reshape((pnm.shape[0], -1))
    expphi = np.exp(1j*mv*phiM)
    Y = pnm*expphi
    return Y