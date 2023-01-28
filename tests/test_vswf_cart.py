import numpy as np
from pytweezer.utils import combined_index, wigner_rotation_matrix
from pytweezer.utils.rotations import rotx, rotz

def test_wigner_rotation_matrix_eye():
    decimal = 3

    val = wigner_rotation_matrix(1, np.eye(3))
    np.testing.assert_array_almost_equal(val, np.eye(3), decimal=decimal, err_msg='Error in the 3x3 eye matrix')
    
    #% Test Nmax = 2
    Nmax = 2
    val = wigner_rotation_matrix(Nmax, np.eye(3))
    num_vals = combined_index(Nmax, Nmax)
    np.testing.assert_array_almost_equal(val, np.eye(num_vals), decimal=decimal, err_msg='Error in the 8x8 eye matrix')


    #% Test Nmax = 20
    Nmax = 20
    val = wigner_rotation_matrix(Nmax, np.eye(3))
    num_vals = combined_index(Nmax, Nmax)
    np.testing.assert_array_almost_equal(val, np.eye(num_vals), decimal=decimal, err_msg='Error in the 8x8 eye matrix')


def test_wigner_rotation_matrix():

    decimal = 3
    
    #  % Results for Nmax=3, R = rotx(32)*rotz(18))
    n1 = np.array([[0.8788 - 0.2855j, 0.1158 + 0.3564j, -0.0723 + 0.0235j],
        [0.0000 + 0.3747j, 0.8480 + 0.0000j, 0.0000 + 0.3747j],
       [-0.0723 - 0.0235j, -0.1158 + 0.3564j, 0.8788 + 0.2855j]]);
    n2 = np.array([[0.6908 - 0.5019j, 0.2878 + 0.3961j, -0.1391 + 0.1011j, -0.0237 - 0.0326j, 0.0047 - 0.0034j],
        [0.1513 + 0.4657j, 0.6117 - 0.1988j, 0.1701 + 0.5235j, -0.1948 + 0.0633j, -0.0124 - 0.0383j],
       [-0.1720 + 0.0000j, 0.0000 + 0.5504j, 0.5788 + 0.0000j, 0.0000 + 0.5504j,  -0.1720 + 0.0000j],
        [0.0124 - 0.0383j, -0.1948 - 0.0633j, -0.1701 + 0.5235j, 0.6117 + 0.1988j, -0.1513 + 0.4657j],
        [0.0047 + 0.0034j, 0.0237 - 0.0326j, -0.1391 - 0.1011j, -0.2878 + 0.3961j, 0.6908 + 0.5019j]])
    n3 = np.array([[0.4637 - 0.6383j, 0.4483 + 0.3257j, -0.1477 + 0.2033j,  -0.0673 - 0.0489j, 0.0121 - 0.0167j, 0.0030 + 0.0022j, -0.0003 + 0.0004j],
        [0.3257 + 0.4483j, 0.3759 - 0.2731j,   0.3513 + 0.4836j, -0.2638 + 0.1917j, -0.0663 - 0.0913j, 0.0212 - 0.0154j, 0.0022 + 0.0030j],
       [-0.2389 + 0.0776j,   0.1847 + 0.5685j,   0.2872 - 0.0933j, 0.1841 + 0.5665j, -0.3300 + 0.1072j, -0.0349 - 0.1073j, 0.0196 - 0.0064j],
        [0.0000 - 0.0832j, -0.3261 + 0.0000j,   0.0000 + 0.5957j, 0.2527 + 0.0000j, 0.0000 + 0.5957j, -0.3261 + 0.0000j, 0.0000 - 0.0832j],
        [0.0196 + 0.0064j,   0.0349 - 0.1073j,  -0.3300 - 0.1072j, -0.1841 + 0.5665j, 0.2872 + 0.0933j, -0.1847 + 0.5685j, -0.2389 - 0.0776j],
       [-0.0022 + 0.0030j,   0.0212 + 0.0154j,   0.0663 - 0.0913j, -0.2638 - 0.1917j, -0.3513 + 0.4836j, 0.3759 + 0.2731j, -0.3257 + 0.4483j],
       [-0.0003 - 0.0004j,  -0.0030 + 0.0022j,   0.0121 + 0.0167j, 0.0673 - 0.0489j, -0.1477 - 0.2033j, -0.4483 + 0.3257j, 0.4637 + 0.6383j]])

    Nmax = 3
    val = wigner_rotation_matrix(Nmax, np.matmul(rotx(32), rotz(18)))
    np.testing.assert_array_almost_equal(val[:3, :3], n1, decimal=decimal, err_msg='Error in the 3x3 matrix')
    np.testing.assert_array_almost_equal(val[3:8, 3:8], n2, decimal=decimal, err_msg='Error in the 5x5 matrix')
    np.testing.assert_array_almost_equal(val[8:, 8:], n3, decimal=decimal, err_msg='Error in the final block matrix')