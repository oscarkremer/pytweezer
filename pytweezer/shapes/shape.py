
from abc import ABC, abstractmethod
import numpy as np
import pickle
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
    if name.lower() not in ('sphere', 'cylinder', 'ellipsoid', 
    'superellipsoid', 'cone-tipped-cylinder', 'cube', 'axisym', 'obj', 
    'stl'):
      raise ValueError('Unsupported simple particle shape')
    elif name == 'sphere':
      shape = Sphere(**kwargs)
    elif name == 'cylinder':
      shape = Cylinder(**kwargs)
    elif name == 'ellipsoid':
      shape = Ellipsoid(**kwargs)
    elif name == 'superellipsoid':
      shape = Superellipsoid(**kwargs)
    elif name == 'cone-tipped-cylinder':
      radius = kwargs['radius']
      height = kwargs['height']
      cone_height = kwargs[2]
      z = np.array([height/2 + cone_height, height/2, 
      -height/2, -height/2 - cone_height])
      rho = np.array([ 0.0, radius, radius, 0.0 ])
      shape = AxiSymLerp(rho, z)
    elif name == 'cube':
      shape = Cube(parameters(:));
    elif name == 'axisym':
      shape = ott.shapes.AxiSymLerp(**kwargs)
    elif name == 'obj'
      shape = WavefrontObj(**kwargs)
    elif name ==  'stl'
      shape = StlLoader(parameters(:));

    return shape    


  @abstractmethod
  def inside(self, shape, varargin):
    pass

  @abstractmethod
  def get_maxRadius(self, shape, varargin):
    pass

  @abstractmethod
  def get_volume(shape, varargin):
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
  

    def voxels(self, shape, varargin):
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
      p = inputParser;
      p.addOptional('spacing', shape.maxRadius/10);
      p.addParameter('plotoptions', []);
      p.addParameter('visualise', nargout == 0);
      p.addParameter('scale', 1.0);
      p.addParameter('axes', []);
      p.addParameter('origin', 'shape');
      p.addParameter('even_range', false);
      p.parse(varargin{:});

      plotoptions = p.Results.plotoptions;
      if isempty(plotoptions)
        plotoptions = {...
          'MarkerFaceColor', 'w', ...
          'MarkerEdgeColor', [.5 .5 .5], ...
          'MarkerSize', 20*p.Results.spacing/shape.maxRadius/p.Results.scale};
      end

      % Calculate range of dipoles
      numr = ceil(shape.maxRadius * p.Results.scale / p.Results.spacing);

      % Add an extra point so we don't have a point around zero
      if p.Results.even_range
        numr = numr + 0.5;
      end

      rrange = (-numr:numr)*p.Results.spacing;

      % Generate the voxel grid
      [xx, yy, zz] = meshgrid(rrange, rrange, rrange);

      % Determine which points are inside
      mask = shape.insideXyz(xx / p.Results.scale, ...
          yy / p.Results.scale, zz / p.Results.scale, ...
          'origin', 'shape');
      xyz = [xx(mask).'; yy(mask).'; zz(mask).'];

      % Translate to world origin
      if strcmpi(p.Results.origin, 'world')
        xyz = xyz + shape.position;
      elseif strcmpi(p.Results.origin, 'shape')
        % Nothing to do
      else
        error('origin must be ''world'' or ''shape''''''); 
      end

      % Visualise the result
      if p.Results.visualise
        
        % Get the axes to use
        our_axes = p.Results.axes;
        if isempty(our_axes)
          our_axes = axes();
        end
        
        plot3(our_axes, xyz(1,:), xyz(2,:), xyz(3,:), 'o', plotoptions{:});
        axis(our_axes, 'equal');
        title(our_axes, ['spacing = ' num2str(p.Results.spacing) ...
            ', N = ' int2str(sum(mask))])
      end

      % Assign output
      if nargout ~= 0
        varargout = {xyz};
      end
    return varargout

  '''
  def insideXyz(self, shape, x, varargin)
      % INSIDEXYZ determine if Cartesian point is inside the shape
      %
      % b = inside(shape, x, y, z) determine if the Cartesian point
      % [x, y, z] is inside the star shaped object.
      %
      % b = insideXyz(shape, xyz) as above, but using a 3xN matrix of
      % [x; y; z] positions.
      %
      % Optional arguments
      %   - origin (enum) -- Coordinate system origin.  Either 'world'
      %     or 'shape' for world coordinates or shape coordinates.
      %
      % See also INSIDE.

      p = inputParser;
      p.addOptional('y', [], @isnumeric);
      p.addOptional('z', [], @isnumeric);
      p.addParameter('origin', 'world');
      p.parse(varargin{:});
      
      if isempty(p.Results.y) && isempty(p.Results.z)
        y = x(2, :);
        z = x(3, :);
        x = x(1, :);
      elseif ~isempty(p.Results.y) && ~isempty(p.Results.z)
        x = x(:);
        y = p.Results.y(:);
        z = p.Results.z(:);
        [x, y, z] = ott.utils.matchsize(x, y, z);
      else
        error('Must suply either 3xN matrix or x, y and z');
      end

      % Translate to shape origin
      if strcmpi(p.Results.origin, 'world')
        x = x - shape.position(1);
        y = y - shape.position(2);
        z = z - shape.position(3);
      elseif strcmpi(p.Results.origin, 'shape')
        % Nothing to do
      else
        error('origin must be ''world'' or ''shape'''''');
      end

      % Convert to spherical coordinates
      [r, t, p] = ott.utils.xyz2rtp(x, y, z);

      % Call the spherical coordinate version
      b = shape.inside(r, t, p, 'origin', 'shape');
    end

    function shape = set.position(shape, value)
      % Check position values
      assert(numel(value) == 3, 'Position must be 3 element vector');
      assert(isnumeric(value), 'Position must be numeric');
      shape.position = value(:);
    end
  end