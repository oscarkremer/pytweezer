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
    J = np.array([[np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)],
            [np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi),-np.sin(theta)],
            [-np.sin(phi), np.cos(phi), np.zeros(theta.shape[0])]])
    print(J)
    xyzv=[xv,yv,zv]
    rv = dot(J[1:length(theta),:],xyzv,2)
    thetav = dot(J[length(theta)+1:2*length(theta),:],xyzv,2);
    phiv = dot(J[2*length(theta)+1:3*length(theta),:],xyzv,2);
    return rv, thetav, phiv, r, theta, phi
