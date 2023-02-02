from .shape import Shape
import pickle
import numpy as np
from abc import ABC, abstractmethod
from pytweezer.utils import matchsize, xyz2rtp, rtp2xyz
from .cylinder import Cylinder
from .ellipsoid import Ellipsoid
from .sphere import Sphere
from .superellipsoid import SuperEllipsoid
from .shape import Shape


class StarShape(Shape, ABC):
    def __init__(self):
        pass

    @abstractmethod
	def radii(self, shape, theta, phi):
		pass

    @abstractmethod
	def normals(self, shape, theta, phi):
		pass

	@abstractmethod
	def axial_symmetry(self, shape):
		pass

    def mirror_symmetry(self, shape):

        axialSym = shape.axialSymmetry()
        orthSym = mod(axialSym, 2) == 0
        mirrorSym = [ orthSym(2) | orthSym(3),
            orthSym(1) | orthSym(3), orthSym(1) | orthSym(2) ]
        return varargout

    def locations(shape, theta, phi):
        theta = theta[:]
        phi = phi[:]
        theta, phi = matchsize(theta, phi)
        varargout{1:nargout} = rtp2xyz(
            shape.radii(theta, phi), theta, phi)
        return varargout
        
    def surf(self, shape, points=[], n_points=[100,100], surf_options=[], position=[], rotation=[], axes=[]):
        if not points:
            sz = npoints
            if sz.size == 1
                sz = np.array([sz sz])
    
            theta, phi = shape.angulargrid('full', true, 'size', sz)
        else:
            theta = points[0]
            phi = points[1]

            if min(theta.shape)==1 and min(phi.shape) == 1:
                phi, theta = np.meshgrid(phi, theta)
            elif size(theta) != size(phi):
                raise ValueError('theta and phi must be vectors or matricies of the same size');    
            sz = theta.shape
    
        X, Y, Z = shape.locations(theta, phi)
        if ~isempty(p.Results.rotation):
            XYZ = [X, Y, Z].'
            XYZ = p.Results.rotation * XYZ;
            X = XYZ(1, :)
            Y = XYZ(2, :)
            Z = XYZ(3, :)
        if ~isempty(p.Results.position)
            X = X + p.Results.position(1)
            Y = Y + p.Results.position(2)
            Z = Z + p.Results.position(3)
        X = reshape(X, sz)
        Y = reshape(Y, sz)
        Z = reshape(Z, sz)
        X(:, end+1) = X(:, 1)
        Y(:, end+1) = Y(:, 1)
        Z(:, end+1) = Z(:, 1)
        if nargout == 0 || ~isempty(p.Results.axes):        
            our_axes = p.Results.axes;
            if isempty(our_axes):
               our_axes = gca();
            
            surf(our_axes, X, Y, Z, p.Results.surfoptions{:});
        varargout = { X, Y, Z }
        return varargout

    def voxels(self, shape, spacing, varargin):
        p = inputParser;
        p.addParameter('plotoptions', {...
            'MarkerFaceColor', 'w', ...
            'MarkerEdgeColor', [.5 .5 .5], ...
            'MarkerSize', 20*spacing/shape.maxRadius});
        p.addParameter('visualise', nargout == 0);
        p.addParameter('origin', 'shape');
        p.addParameter('even_range', false);
        p.parse(varargin{:});
        numr = ceil(shape.maxRadius / spacing);
      
        if p.Results.even_range:
            numr = numr + 0.5
        rrange = (-numr:numr)*spacing;
        [xx, yy, zz] = meshgrid(rrange, rrange, rrange);
        mask = shape.insideXyz(xx, yy, zz, 'origin', 'shape');
        xyz = [xx(mask).T; yy(mask).T; zz(mask).T];
        if strcmpi(p.Results.origin, 'world')
            xyz = xyz + shape.position
        elif strcmpi(p.Results.origin, 'shape')
            pass
        else:
            raise ValueError('Origin must be ''world'' or ''shape''')
        return xyz

    def normalsXyz(self, shape, theta, phi):
        [theta,phi] = ott.utils.matchsize(theta, phi)
        n = rtpv2xyzv(shape.normals(theta, phi), ...
            [ zeros(size(theta)), theta, phi ]);
        return n

    def inside(self, shape, radius, theta, phi, varargin)
      p = inputParser
      p.addParameter('origin', 'world');
      p.parse(varargin{:});

      theta = theta(:);
      phi = phi(:);
      radius = radius(:);
      [radius,theta,phi] = ott.utils.matchsize(radius,theta,phi);

      % Translate to shape origin
      if strcmpi(p.Results.origin, 'world')

        % Only do work if we need to
        if vecnorm(shape.position) ~= 0
          [x,y,z] = ott.utils.rtp2xyz(radius, theta, phi);
          x = x - shape.position(1);
          y = y - shape.position(2);
          z = z - shape.position(3);
          [radius, theta, phi] = ott.utils.xyz2rtp(x, y, z);
        end
      elseif strcmpi(p.Results.origin, 'shape')
        % Nothing to do
      else
        error('origin must be ''world'' or ''shape''');
      end

      assert(all(radius >= 0), 'Radii must be positive')

      % Determine if points are less than shape radii
      r = shape.radii(theta, phi);
      b = radius < r;

    end

    def insideXyz(self, shape, varargin):
        p = inputParser;
        p.addOptional('x', []);
        p.addOptional('y', []);
        p.addOptional('z', []);
        p.addParameter('origin', 'world');
        p.parse(varargin{:});
        if isempty(p.Results.y) && isempty(p.Results.z)
            x = p.Results.x(1, :);
            y = p.Results.x(2, :);
            z = p.Results.x(3, :);
        else
            x = p.Results.x(:);
            y = p.Results.y(:);
            z = p.Results.z(:);
            [x, y, z] = ott.utils.matchsize(x, y, z);
       if strcmpi(p.Results.origin, 'world'):
            x = x - shape.position(1)
            y = y - shape.position(2)
            z = z - shape.position(3)
        elif strcmpi(p.Results.origin, 'shape'):
            pass
        else
            raise ValueError('origin must be ''world'' or ''shape''');
        [r, t, p] = xyz2rtp(x, y, z)
        b = shape.inside(r, t, p, 'origin', 'shape')
        return b

    def angulargrid(self, shape, varargin):

        p = inputParser;
        p.addOptional('Nmax', 100);
        p.addParameter('full', false);    % Not used
        p.addParameter('size', []);    % Not used
        p.parse(varargin{:});

        if isempty(p.Results.size):
            ntheta = 2*(p.Results.Nmax + 2);
            nphi = 3*(p.Results.Nmax + 2) + 1;
            if ~p.Results.full:
                [~, ~, z_axial_symmetry] = shape.axialSymmetry();
                if z_axial_symmetry == 0:
                    ntheta = 4*(p.Results.Nmax + 2);
                    nphi = 1;
                else:
                    nphi = round(nphi / z_axial_symmetry);

                [~, ~, z_mirror_symmetry] = shape.mirrorSymmetry();
                if z_mirror_symmetry:
                    ntheta = round(ntheta / 2)
        else:
            ntheta = p.Results.size(1)
            nphi = p.Results.size(2)
        [~, ~, z_axial_symmetry] = shape.axialSymmetry()
        if ~p.Results.full && z_axial_symmetry == 0:
            nphi = 1
        [theta, phi] = ott.utils.angulargrid(ntheta, nphi)
        if ~p.Results.full
            [~, ~, z_mirror_symmetry] = shape.mirrorSymmetry()
            if z_mirror_symmetry:
                theta = theta / 2.0
            if z_axial_symmetry > 1:
                phi = phi / z_axial_symmetry
        r = shape.radii(theta, phi)
        return r, theta, phi
