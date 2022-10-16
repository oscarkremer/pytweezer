import pytest
import numpy as np
from scipy.sparse import csr_matrix
from pytweezer.utils import combined_index, translate_z

def test_translate_z_eye_values():   
    decimal = 3    
    A, B, C = translate_z(5, 0.0)
    sz = combined_index(5, 5)
    assert A.shape == (sz, sz), 'Error in matrix A dimensions!'
    np.testing.assert_array_almost_equal(A, np.eye(sz, sz), decimal=decimal, 
        err_msg='Error in the A matrix for eye case!')
    np.testing.assert_array_almost_equal(B, np.zeros((sz, sz)), decimal=decimal,
        err_msg='Error in the B matrix for eye case!')
    np.testing.assert_array_almost_equal(C, np.eye(sz, sz), decimal=decimal,
        err_msg='Error in the C matrix for eye case!')
    A, B, C = translate_z(5, np.pi/2)
    targetA = csr_matrix((np.array([-0.0786, 0.1411, 0.1964, -0.0217, -0.2398]), 
            (np.array([0, 0, 0, 0, 0]), np.array([0, 4, 10, 18, 28])), 
            ), shape=(1, sz)).toarray()
    targetB = csr_matrix((np.array([-0.1306j, -0.1357j, 0.1181j, 0.2770j, 0.1312j]),
            (np.array([0, 0, 0, 0, 0]), np.array([0, 4, 10, 18, 28]))), 
            shape=(1, sz)).toarray()
    np.testing.assert_array_almost_equal(A[0,:].reshape((1, -1)), targetA, decimal=decimal, err_msg='Error for A matrix')
    np.testing.assert_array_almost_equal(B[0,:].reshape((1, -1)), targetB, decimal=decimal, err_msg='Error for B matrix')

def test_translate_z_for_n_array():
    decimal = 3
    Afull, Bfull, Cfull = translate_z(7, np.pi/2)
    sz1 = combined_index(5, 5)
    sz2 = combined_index(7, 7)
    A, B, C = translate_z(np.array([5, 7]), np.pi/2)
    assert A.shape == (sz1, sz2), 'Incorrect matrix size [5 7]'
    np.testing.assert_array_almost_equal(A, Afull[:sz1, :sz2], decimal=decimal, err_msg='Error for A matrix')
    np.testing.assert_array_almost_equal(B, Bfull[:sz1, :sz2], decimal=decimal, err_msg='Error for B matrix')
    A, B, C = translate_z(np.array([7, 5]), np.pi/2)
    assert A.shape == (sz2, sz1), 'Incorrect matrix size [7 5]'
    np.testing.assert_array_almost_equal(A, Afull[:sz2, :sz1], decimal=decimal, err_msg='Error for A matrix')
    np.testing.assert_array_almost_equal(B, Bfull[:sz2, :sz1], decimal=decimal, err_msg='Error for B matrix')



'''
function testLargeTranslations(testCase)

  import matlab.unittest.constraints.IsEqualTo;
  import matlab.unittest.constraints.AbsoluteTolerance;
  tol = 1.0e-6;

  % For this translation the beam power should go to zero

  beam = ott.BscPmGauss();
  beam.power = 1.0;
  tbeam = beam.translateXyz([300;0;0]);  % calls translate_z

  testCase.verifyThat(tbeam.power, IsEqualTo(0.0, ...
      'Within', AbsoluteTolerance(tol)), ...
      'Beam power does not drop to zero for large radial translations');

end
'''
def test_negative_translate_z():
    Nmax = 1
    decimal = 3 
    A, B, C = translate_z(Nmax, 1e-6, method='gumerov')
    expected_A = np.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
    expected_B = np.array([[0.0+3.14159265e-06j, 0, 0],
                [0, 0, 0],
                [0, 0, -3.14159265e-06j]]) 
    expected_C = np.array([[[1.0,  0.0],
                            [0.0, 0.0]],
                            [[0.0,  0.0],
                            [ 1.0, 1.0]]])
    np.testing.assert_array_almost_equal(A, expected_A, decimal=decimal, err_msg='Gumerov failed for negative translation')
    np.testing.assert_array_almost_equal(B, expected_B, decimal=decimal, err_msg='Gumerov failed for negative translation')
    np.testing.assert_array_almost_equal(C, expected_C, decimal=decimal, err_msg='Gumerov failed for negative translation')  
    A, B, C = translate_z(Nmax, 1e-6, method='videen')
    expected_A = np.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
    expected_B = np.array([[0.0+3.14159265e-06j, 0, 0],
                [0, 0, 0],
                [0, 0, -3.14159265e-06j]]) 
    expected_C = np.array([[[1.0,  0.0],
                            [0.0, 0.0]],
                            [[0.0,  0.0],
                            [ 1.0, 1.0]]])
    np.testing.assert_array_almost_equal(A, expected_A, decimal=decimal, err_msg='Videen failed for negative translation')
    np.testing.assert_array_almost_equal(B, expected_B, decimal=decimal, err_msg='Videen failed for negative translation')
    np.testing.assert_array_almost_equal(C, expected_C, decimal=decimal, err_msg='Videen failed for negative translation')

 
def test_negative_translation():
    Nmax = 1
    decimal = 3 
    A, B, C = translate_z(Nmax, -np.pi/2, method='gumerov', function_type='sbesselh1')
    expected_A = np.array([[-0.0393135 +0.06457856j, 0, 0],[0, 0.01322905+0.0080346j, 0],[0, 0,  -0.0393135 +0.06457856j]])

    expected_B = np.array([[-0.03964919+0.06528276j , 0, 0],
                [0, 0, 0],
                [0, 0, 0.03964919-0.06528276j]]) 
    expected_C = np.array([[[-0.02179931+0.04573057j, 0+0.j],
                            [-0.07538204-0.04578294j,0+0.j]],
                            [[0.07538204+0.04578294j, 0+0.j],
                            [-0.09185605+0.12112251j, 0.01322905+0.0080346j]]])
    np.testing.assert_array_almost_equal(A, expected_A, decimal=decimal, err_msg='Gumerov failed for negative translation')
    np.testing.assert_array_almost_equal(B, expected_B, decimal=decimal, err_msg='Gumerov failed for negative translation')
    np.testing.assert_array_almost_equal(C, expected_C, decimal=decimal, err_msg='Gumerov failed for negative translation')  
    A, B, C = translate_z(Nmax, -np.pi/2, method='videen', function_type='sbesselh2')
 
    expected_A = np.array([[-0.078627-0.12915712j, 0, 0],[0,  0.02645811-0.01606921j, 0],[0, 0, -0.078627  -0.12915712j]])
    expected_B = np.array([[ 0.07929837+0.13056553j, 0, 0],
                [0, 0, 0],
                [0, 0,  -0.07929837-0.13056553j]]) 
    expected_C = np.array([[[-0.04359863-0.09146115j, 0+0j],
                            [-0.15076408+0.09156587j, 0+0j]],
                        [[ 0.15076408-0.09156587j,0+0j],
                        [-0.1837121-0.24224503j,0.02645811-0.01606921j]]])
    np.testing.assert_array_almost_equal(A, expected_A, decimal=decimal, err_msg='Videen failed for negative translation')
    np.testing.assert_array_almost_equal(B, expected_B, decimal=decimal, err_msg='Videen failed for negative translation')
    np.testing.assert_array_almost_equal(C, expected_C, decimal=decimal, err_msg='Videen failed for negative translation')

def test_wrong_method_parameter():
    Nmax = 1
    with pytest.raises(ValueError):
        A, B, C = translate_z(Nmax, -np.pi/2, method='random', function_type='sbesselh1')
    
def test_wrong_function_type_parameter():
    Nmax = 1
    with pytest.raises(ValueError):
        A, B, C = translate_z(Nmax, -np.pi/2, method='random', function_type='legendre')
