import pytest
import numpy as np
from pytweezer.beams import Gaussian

def test_sphere_perimeter():
    beam = Gaussian()
    beam.power = 1.0
    dz = np.pi/2

    translated_beam1 = beam.translate_z(z=dz)

    _, AB = beam.translate_z(dz)
#  tbeam2 = AB * beam;

#  [~, A, B] = beam.translateZ(dz);
#  tbeam3 = [ A B ; B A ] * beam;
#  tbeam4 = beam.translate(A, B);

#  % Check all beams are equal
#  target = tbeam1.getCoefficients;
#  testCase.verifyThat(tbeam2.getCoefficients, IsEqualTo(target, ...
#    'Within', AbsoluteTolerance(tol)), ...
#    'AB * beam does not match translateZ');
#  testCase.verifyThat(tbeam3.getCoefficients, IsEqualTo(target, ...
#    'Within', AbsoluteTolerance(tol)), ...
#    '[A, B; B, A] * beam does not match translateZ');
#  testCase.verifyThat(tbeam4.getCoefficients, IsEqualTo(target, ...
#    'Within', AbsoluteTolerance(tol)), ...
#    'beam.translate(A, B) does not match translateZ');
#
#  xbeam1 = beam.translateXyz([dz; dz; 0]);
#
#  [~, AB] = beam.translateXyz([dz; dz; 0]);
#  xbeam2 = AB * beam;
#
#  [~, A, B] = beam.translateXyz([dz; dz; 0]);
#  xbeam3 = beam.translate(A, B);
#
#  [~, A, B, D] = beam.translateXyz([dz; dz; 0]);
#  xbeam4 = beam.translateXyz(A, B, D);
#
#  % Check all beams are equal
#  target = xbeam1.getCoefficients;
#  testCase.verifyThat(xbeam2.getCoefficients, IsEqualTo(target, ...
#    'Within', AbsoluteTolerance(tol)), ...
#    'AB * beam does not match translateXyz');
#  testCase.verifyThat(xbeam3.getCoefficients, IsEqualTo(target, ...
#    'Within', AbsoluteTolerance(tol)), ...
#    'beam.translate(A, B) does not match translateXyz');
#  testCase.verifyThat(xbeam4.getCoefficients, IsEqualTo(target, ...
#    'Within', AbsoluteTolerance(tol)), ...
#    'beam.translateXyz(A, B, D) does not match translateXyz');


'''
    decimal = 3
    radius = 1 
    shape = Sphere(radius)
    np.testing.assert_array_almost_equal(shape.perimeter, 6.2832, decimal=decimal, err_msg='Error computing perimeter of sphere')

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
function testUnevenTranslation(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;

  beam = testCase.TestData.beam;
  beam.power = 1.0;
  dz = pi/2;
  tol = 1.0e-6;

  tbeam = beam.translateZ(dz, 'Nmax', beam.Nmax - 5);
  testCase.verifyThat(tbeam.Nmax, IsEqualTo(beam.Nmax - 5), ...
    'Translated beam does not have correct Nmax');

end

function testMakeBeamVectorEmpty(testCase)
% Check that make_beam_vector works for empty beams

  import matlab.unittest.constraints.IsEqualTo;

  a = [];
  b = [];
  nn = [];
  mm = [];
  
  [a1, b1] = ott.Bsc.make_beam_vector(a, b, nn, mm);
  
  testCase.verifyThat([size(a1), size(b1)], IsEqualTo([0, 0, 0, 0]), ...
    'Size of beam vectors incorrect with implicit Nmax');
  
  [a2, b2] = ott.Bsc.make_beam_vector(a, b, nn, mm, 1);
  
  testCase.verifyThat([size(a2), size(b2)], IsEqualTo([3, 0, 3, 0]), ...
    'Size of beam vectors incorrect with explicit Nmax');
end

function testMakeBeamVectorMulti(testCase)
% Check to make sure make_beam_vector functions with multiple beams
% with the same nn and mm indices.

  a = [1; 2; 3];
  b = [4; 5; 6];
  
  nn = [1; 2; 3];
  mm = [0; 0; 0];
  
  [a1, b1] = ott.Bsc.make_beam_vector(a, b, nn, mm);
  [a2, b2] = ott.Bsc.make_beam_vector(a+6, b+6, nn, mm);
  
  [ac, bc] = ott.Bsc.make_beam_vector([a, a+6], [b, b+6], nn, mm);
  
  import matlab.unittest.constraints.IsEqualTo;
  
  testCase.verifyThat([ac(:, 1); bc(:, 1)], IsEqualTo([a1; b1]), ...
    'First beam doesn''t match');
  
  testCase.verifyThat([ac(:, 2); bc(:, 2)], IsEqualTo([a2; b2]), ...
    'Second beam doesn''t match');

end

function testSum(testCase)

  beam1 = ott.BscPmGauss('polarisation', [0, 1i]);
  beam2 = ott.BscPmGauss('polarisation', [1, 0]);
  beam3 = beam1 + beam2;
  
  beamU = beam1.append(beam2);
  testCase.verifyEqual(beamU.sum(), beam3, ...
    'beamU.sum() incorrect');
  testCase.verifyEqual(sum(beamU), beam3, ...
    'sum(beamU) incorrect');
  
  % Test array sum
  beamarr = [beam1, beam2];
  testCase.verifyEqual(sum(beamarr), beam3, ...
    'sum([beam1, beam2]) incorrect');
  
end
'''
