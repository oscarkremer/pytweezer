import pickle
import numpy as np
from abc import ABC, abstractmethod
from pytweezer.utils import match_size, xyz2rtp

class Shape(ABC):
    '''
    %Shape abstract class for optical tweezers toolbox shapes
    %
    % Properties
    %   maxRadius     maximum distance from shape origin
    %   volume        volume of shape
    %   position      Location of shape ``[x, y, z]``
    %
    % Methods (abstract):
    %   inside(shape, ...) determine if spherical point is inside shape
    %
    % Methods:
    %   writeWavefrontObj(shape, ...) write shape to Wavefront OBJ file
    %       only implemented if shape supports this action.
    %   insideXyz(shape, ...) determine if Cartesian point is inside shape
    %       requires inside(shape, ...) to be implemented.
    %   simple(...) simplified constructor for shape-like objects.
    %
    % See also simple, ott.shapes.Cube, ott.shapes.TriangularMesh.

    % This file is part of the optical tweezers toolbox.
    See LICENSE.md for information about using/distributing this file.
    '''

    def __init__(self, position: np.array = np.array([0, 0, 0])):
        '''
            properties
            position = [0;0;0];       % Location of shape ``[x, y, z]``
        '''
        self.position = position

    @abstractmethod
    def inside(self, shape):
        pass

    @abstractmethod
    def get_max_radius(self, shape):
        pass

    @abstractmethod
    def get_volume(self, shape):
        pass

    def _save_(self, filename, verts, faces):
        '''% Helper to write faces and vertices to a file
        %
        % helper(filename, verts, faces)
        %   verts 3xN array of floats
        %   faces mxN array of integers (starting at 1)
        '''
        with open(f'verts_{filename}.pickle', 'wb') as verts_file:
            pickle.dump(verts, verts_file)
        with open(f'faces_{filename}.pickle', 'wb') as faces_file:
            pickle.dump(faces, faces_file)                     
        with open('filename.pickle', 'rb') as handle:
            b = pickle.load(handle)
      
    def writeWavefrontObj(self, filename):
        '''
        Write representation of shape to Wavefront OBJ file
        %
        % This is a placeholder, not supported by all shape types
        '''
        raise TypeError('Shape type does not support writing to WavefrontOBJ file')

    def surf(self, varargin):
        '''
        % SURF generate a visualisation of the shape
        %
        % SURF() displays a visualisation of the shape in the current figure.
        %
        % SURF(..., 'surfoptions', {varargin}) specifies the options to
        % pass to the surf function.
        '''
        raise TypeError('Shape type does not support surf visualisation');
	

    def voxels(self, shape, spacing=0, plotoptions=[], visualise=False, plot=False, 
            scale=1.0, axes=[], origin='shape', even_range=False, inline=True):
        '''
        % Generate an array of xyz coordinates for voxels inside the shape
        %
        % Usage
        %   voxels(spacing) shows a visualisation of the shape with
        %   circles placed at locations on a Cartesian grid.
        %
        %   xyz = voxels(spacing) returns the voxel locations.
        %
        % Optional named arguments
        %   - 'plotoptions'   Options to pass to the plot3 function
        %   - 'visualise'     Show the visualisation (default: nargout == 0)
        %   - origin (enum) -- Coordinate system origin.  Either 'world'
        %     or 'shape' for world coordinates or shape coordinates.
        '''
        if not spacing:
            spacing = shape.maxRadius/10
        if not plotoptions:
            plotoptions = {'markerfacecolor': 'w',
                        'markeredgecolor': [0.5, 0.5, 0.5], 
                        'markersize': 20*spacing/shape.maxRadius/Results.scale}
        if not inline:
            plot=True
        numr = np.ceil(shape.maxRadius * scale / spacing)
        if even_range:
            numr = np.arange() + 0.5
        rrange = np.arange(-numr, numr+1)*spacing
        xx, yy, zz = np.meshgrid(rrange, rrange, rrange)
        mask = shape.insideXyz(xx/scale, yy/scale, zz/scale, 'origin', 'shape')
        xyz = [xx(mask).T, yy(mask).T, zz(mask).T]
        if strcmpi(origin, 'world'):
            xyz = xyz + shape.position
        elif strcmpi(origin, 'shape'):
            pass
        else:
            raise ValueError('origin must be `world` or `shape`')
        if plot:
            pass
        if not inline:
            return insideXyz

    def inside_xyz(self, shape, x, y=np.array([]), z=np.array([]), origin='world'):
        if not y or not z:
            y = x[1, :]
            z = x[2, :]
            x = x[1, :]
        elif y.shape and z.shape:
            x = x.T
            y = z.T
            z = z.T
            x, y, z = match_size(x, y, z)
        else:
            raise ValueError('Must suply either 3xN matrix or x, y and z')
        if origin.lower() == 'world':
            x = x - shape.position[0]
            y = y - shape.position[1]
            z = z - shape.position[2]
        elif origin.lower() == 'shape':
            pass
        else:
            raise ValueError('Origin must be `world` or `shape`')
        r, t, p = xyz2rtp(x, y, z)
        return shape.inside(r, t, p, origin='shape')