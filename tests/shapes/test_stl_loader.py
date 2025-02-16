
import pytest
import numpy as np
from pytweezer.shapes import STLLoader

def testLoad():
	shape = STLLoader('cube.stl')
	print(shape.faces, shape.verts)
#	  testCase.verifyThat(numel(shape.verts), IsEqualTo(3*8), ...
#    'Incorrect number of vertices');

#  testCase.verifyThat(numel(shape.faces), IsEqualTo(3*2*6), ...
#    'Incorrect number of faces');

'''
def test_sphere_volume():
    decimal = 3
    radius = 1 
    shape = Sphere(radius)
    np.testing.assert_array_almost_equal(shape.volume, 4.1888, decimal=decimal, err_msg='Error computing volume of sphere')


def test_sphere_inside():
    radius = 1
    shape = Sphere(radius)
    x = np.append(0.5*1*np.random.rand(1, 3), 4.0)
    y = np.append(0.5*1*np.random.rand(1, 3), 4.0)
    z = np.append(0.5*1*np.random.rand(1, 3), 4.0)
    for value, expected in zip(shape.inside_xyz(x,y,z), [True, True, True, False]):
        assert value==expected, 'Error in inside_xyz method for Spherical Shape'

'''

'''
function tests = stlloader
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)
  addpath('../../');
end

function testLoad(testCase)

  % Load the shape
  shape = ott.shapes.StlLoader('cube.stl');

  import matlab.unittest.constraints.IsEqualTo;

  testCase.verifyThat(numel(shape.verts), IsEqualTo(3*8), ...
    'Incorrect number of vertices');

  testCase.verifyThat(numel(shape.faces), IsEqualTo(3*2*6), ...
    'Incorrect number of faces');

end

function testInsideXyz(testCase)

  shape = ott.shapes.StlLoader('cube.stl');
  radius = shape.maxRadius;
  
  % Choose three points inside the shape and one outside
  b = [true, true, true, false].';
  x = [0.5*radius.*rand(1, 3), 4.0];
  y = [0.5*radius.*rand(1, 3), 4.0];
  z = [0.5*radius.*rand(1, 3), 4.0];
  
  xyz = [x(:), y(:), z(:)].';
  
  testCase.verifyEqual(shape.insideXyz(x, y, z), b, ...
    'insideXyz with 3 arguments failed');
  testCase.verifyEqual(shape.insideXyz(xyz), b, ...
    'insideXyz with 1 argument failed');
  
  testCase.verifyEqual(shape.insideXyz(x, y, z, 'origin', 'world'), b, ...
    'insideXyz with 3 arguments failed and optional arg');
  testCase.verifyEqual(shape.insideXyz(xyz, 'origin', 'world'), b, ...
    'insideXyz with 1 argument failed and optional arg');
end
'''
