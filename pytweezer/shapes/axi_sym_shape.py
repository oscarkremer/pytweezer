import pickle
import numpy as np
from abc import ABC, abstractmethod
from pytweezer.utils import match_size, xyz2rtp
from .shape import Shape


class AxiSymShape(Shape, ABC):    
    def __init__(self, perimeter):
        self.perimeter = perimeter

    @abstractmethod
    def radii(self, shape, theta, phi):
        pass

    @abstractmethod
    def normals(self, shape, theta, phi):
        pass

    @abstractmethod
    def axial_symmetry(self, shape):
        pass

    @abstractmethod
    def get_perimeter(self, **kwargs):
        pass

    def boundary_points_npts(self, shape, npts=[], n_max=[]):
        if not (npts or n_max):
            raise ValueError('Must specify either npts or Nmax');
        elif npts and n_max:
            raise ValueError('Both number of points and Nmax specified')
        else:
            if npts:
                return npts
            else:
                npts = np.ceil(combined_index(n_max, n_max)**(2/20)+5)
                return npts

    def boundary_points_area(self, shape, rho, z, rho_out, z_out, rtp):
        dst = np.zeros((rtp.shape[0],3))
        dst[1:-1,0] = (rho_out[2:] - rho_out[1:end-1])/2
        dst[1:-1,2] = (z_out[2:end] - z_out[1:end-1])/2
        dst[0, 0] = rho_out[0:2].mean()-rho[0]
        dst[0, 2] = z_out[1:2].mean() - z[0]
        dst[-1, 0] = rho[-1] - rho_out[-2:].mean()
        dst[-1, 2] = z[-1] - z_out[-2:].mean()
        ds = rtp[:,0]*np.sqrt(np.power(np.abs(dst),2).sum(axis=1))*np.sin(rtp[:,1])
        return ds

    def boundary_points_rhoz(self, shape, rho, z, npts=[], n_max=[]):
        ntheta = shape.boundary_points_npts()
        axisym = shape.axial_symmetry()
        if axisym(3) == 0:
            raise ValueError('Only supports axisymetric particles');
        '''    
        ds = shape.perimiter / ntheta / 2.0
        s=sqrt((rho(2:end)-rho(1:end-1)).^2+(z(2:end)-z(1:end-1)).^2);
        zout=zeros(ntheta,1);
        rhoout=zout;
        nxyz=zeros(ntheta,3);
        sdeficit=0;
        ncum=0;
        for ii=2:length(rho)
            N=s(ii-1)/ds;
            Nused=round(N+sdeficit);
            nc=[-(z(ii)-z(ii-1)),0,(rho(ii)-rho(ii-1))];
            nc=nc/norm(nc,2);
            if Nused>=1
                drho=(rho(ii)-rho(ii-1))/N*ones(Nused,1);
                rhot=cumsum(drho)-drho/2-sdeficit*drho(1);
                rhoout(ncum+(1:Nused))=rho(ii-1)+rhot;

                dz=(z(ii)-z(ii-1))/N*ones(Nused,1);

                zt=cumsum(dz)-dz/2-sdeficit*dz(1);
                zout(ncum+(1:Nused))=z(ii-1)+zt;

                nxyz(ncum+(1:Nused),:)=repmat(nc,[length(zt),1]);

                sdeficit=(N-Nused+sdeficit);
            else
                sdeficit=sdeficit+N;
            ncum=ncum+Nused;
        if ncum < ntheta
            warning('OTT:SHAPES:AXISYMSHAPE:boundarypoints_length', ...
            'Number of points generated does not match request');
            zout = zout(1:ncum);
            rhoout = zout(1:ncum);
            nxyz = nxyz(1:ncum, :);
        [n,rtp]=ott.utils.xyzv2rtpv(nxyz,[rhoout,zeros(size(rhoout)),zout]);
        ds = shape.boundarypoints_area(rho, z, rhoout, zout, rtp);
        return [rtp, n, ds]
    '''
    
    def boundary_points(self, shape, varargin):
        pass
        '''
        npts = shape.boundary_points_npts(varargin{:});
        [theta, phi] = ott.utils.angulargrid(npts*2, 1);
        theta = [0.0; theta; pi];
        phi = [phi(1); phi; phi(end)];
        xyz = shape.locations(theta, phi);
        rho = xyz(:, 1);
        z = xyz(:, 3);

        [rtp, n, ds] = shape.boundarypoints_rhoz(rho, z, varargin{:});
        return [rtp, n, ds]
        '''

