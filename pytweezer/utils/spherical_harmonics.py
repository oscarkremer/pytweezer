import numpy as np
from .legendre_row import legendre_row
from .match_size import match_size
from .meshgrid import meshgrid
from .compute_y_exp_minus_and_plus import * 

def spherical_harmonics(n, m, theta, phi):
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
    mi = m
    m = np.arange(-n,n+1)
    index_bigger = np.argwhere(abs(m) <=n)
    if index_bigger.size:
       m = m[tuple(index_bigger.T)]
    theta, phi = match_size(theta, phi)
    pnm = legendre_row(n, theta)
    phiM, mv = meshgrid(phi, m)
    Y, expplus, expminus = compute_y_exp_minus_and_plus_n(m, mv, phi, pnm)
    Y2 = _spherical_harmonics_(n+1,theta,phi) 
    Y, Ytheta, Yphi = compute_y_theta_phi(Y, Y2, expplus, expminus, n, m, mi, mv, theta)
    return Y, Ytheta, Yphi

def _spherical_harmonics_(n, theta, phi):
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
    m = np.arange(-n,n+1)
    index_bigger = np.argwhere(abs(m) <=n)
    if index_bigger.size:
       m = m[tuple(index_bigger.T)]
    theta, phi = match_size(theta, phi)
    pnm = legendre_row(n, theta)
    phiM, mv = meshgrid(phi, m)
    Y = compute_y_exp_n(m, mv, phiM, pnm)
    return Y