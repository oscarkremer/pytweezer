import pickle
import numpy as np
from abc import ABC, abstractmethod
from pytweezer.utils import matchsize, xyz2rtp
from .cylinder import Cylinder
from .ellipsoid import Ellipsoid
from .sphere import Sphere
from .superellipsoid import SuperEllipsoid


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


	@staticmethod
	def simple(name, **kwargs):
		'''
			function shape = simple(name, parameters)
				%SIMPLE constructs a shape for a bunch of simple particles
			% This method replaces shapesurface.m from OTTv1.
			%
			% Supported shapes [parameters]:
			%   'sphere'          Sphere [ radius ]
			%   'cylinder'        z-axis aligned cylinder [ radius height ]
			%   'ellipsoid'       Ellipsoid [ a b c]
			%   'superellipsoid'  Superellipsoid [ a b c e n ]
			%   'cone-tipped-cylinder'      [ radius height cone_height ]
			%   'cube'            Cube [ width ]
			%   'axisym'          Axis-symetric particle [ rho(:) z(:) ]
			%   'stl'             load STL file [filename]
			%   'obj'             load Wavefront OBJ file [filename]
		'''
		name = name.lower()
		if name not in ('sphere', 'cylinder', 'ellipsoid', 
		'superellipsoid', 'cone-tipped-cylinder', 'cube', 'axisym', 'obj', 
		'stl'):
			raise ValueError('Unsupported simple particle shape')
		elif name == 'sphere':
			shape = Sphere(kwargs)
		elif name == 'cylinder':
			shape = Cylinder(kwargs)
		elif name == 'ellipsoid':
			shape = Ellipsoid(kwargs)
		elif name == 'superellipsoid':
			shape = SuperEllipsoid(kwargs)
		elif name == 'cone-tipped-cylinder':
			radius = kwargs['radius']
			height = kwargs['height']
			cone_height = kwargs['cone_height']
			z = np.array([height/2 + cone_height, height/2, 
			-height/2, -height/2 - cone_height])
			rho = np.array([ 0.0, radius, radius, 0.0 ])
			shape = AxiSymLerp(rho, z)
		elif name == 'cube':
			shape = Cube(kwargs)
		elif name == 'axisym':
			shape = AxiSymLerp(kwargs)
		elif name == 'obj'
			shape = WavefrontObj(kwargs)
		elif name ==  'stl'
			shape = StlLoader(kwargs)
		return shape    

	@abstractmethod
	def inside(self, shape):
		pass

	@abstractmethod
	def get_maxRadius(self, shape):
		pass

	@abstractmethod
	def get_volume(shape):
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
      
	def getMaxRadius(self, shape):
		return shape.getMaxRadius()  

	def getVolume(self, shape):
		return shape.getVolume()

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
	

    def voxels(self, shape, spacing=0, plotoptions=[], visualise, plot=False, 
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

      # Generate the voxel grid
		xx, yy, zz = np.meshgrid(rrange, rrange, rrange);

		# Determine which points are inside
		mask = shape.insideXyz(xx/scale, yy/scale, zz/scale,
			'origin', 'shape')
		xyz = [xx(mask).T; yy(mask).T; zz(mask).T]

		# Translate to world origin
		if strcmpi(origin, 'world'):
			xyz = xyz + shape.position
		elif strcmpi(origin, 'shape'):
			#% Nothing to do
			pass
		else:
			raise ValueError('origin must be `world` or `shape`')
		if plot:
			% Get the axes to use
			our_axes = p.Results.axes;
			if isempty(our_axes)
			our_axes = axes();
			end
			
			plot3(our_axes, xyz(1,:), xyz(2,:), xyz(3,:), 'o', plotoptions{:});
			axis(our_axes, 'equal');
			title(our_axes, ['spacing = ' num2str(p.Results.spacing) ...
				', N = ' int2str(sum(mask))])

		% Assign output
		if not inline:
			return insideXyz

	def insideXyz(self, shape, x, y=np.array([]), z=np.array([]), origin='world'):
		if not y or not z:
			y = x[1, :]
			z = x[2, :]
			x = x[1, :]
		elif y.shape and z.shape:
			x = x.T
			y = z.T
			z = z.T
        	x, y, z = matchsize(x, y, z)
	      else:
    		raise ValueError('Must suply either 3xN matrix or x, y and z')
		#		% Translate to shape origin
		if origin.lower() == 'world':
			x = x - shape.position[0]
			y = y - shape.position[1]
			z = z - shape.position[2]
		elif origin.lower() == 'shape':
			pass
		else:
			raise ValueError('Origin must be `world` or `shape`')

#		% Convert to spherical coordinates
		r, t, p = xyz2rtp(x, y, z)

		#% Call the spherical coordinate version
		return shape.inside(r, t, p, origin='shape')