import numpy as np
from .xyz2rtp import xyz2rtp

def xyzv2rtpv(xv, yv, zv, x=np.array([]), y=np.array([]), z=np.array([])):
    '''
    % XYZV2RTPV cartiesian to spherical vector field conversion
    %
    % [rv,thetav,phiv,r,theta,phi] = XYZV2RTPV(xv,yv,zv,x,y,z)
    %
    % [vec_sph,pos_sph] = XYZV2RTPV(vec_cart,pos_cart)
    %
    % See also rtpv2xyzv and xyz2rtp.

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.
    '''
    r, theta, phi = xyz2rtp(x, y, z)
    J1 = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]).T
    J2 = np.array([np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi),-np.sin(theta)]).T
    J3 = np.array([-np.sin(phi), np.cos(phi), np.zeros(theta.shape[0])]).T
    J = np.concatenate([J1, J2, J3], axis=0)
    xyzv = np.array([xv, yv, zv]).T
    rv = np.sum(J[:theta.shape[0],:].conj()*xyzv, axis=1)
    thetav = np.sum(J[theta.shape[0]:2*theta.shape[0],:].conj()*xyzv, axis=1)
    phiv = np.sum(J[2*theta.shape[0]:,:].conj()*xyzv, axis=1)
    return rv, thetav, phiv, r, theta, phi
