
from pytweezer.utils.rotations import *
from pytweezer.utils import combined_index, wigner_rotation_matrix
from copy import copy
def rotate_yz(beam, angley, anglez, r, n_max=0):
    n_max = beam.get_n_max() if not n_max else n_max
    R = np.matmul(rotz(anglez*180/np.pi),roty(angley*180/np.pi))
    beam, D = rotate(beam, R=R)
    return beam, D


def rotate(beam, wigner=np.array([]), R=np.array([]), n_max=0):
    n_max = beam.get_n_max() if not n_max else n_max
    if R.size and not wigner.size:        
        if np.power(np.eye(3)-R, 2).sum() < 1e-6:
            D = np.eye(combined_index(n_max, n_max))
            return beam, D
        D = wigner_rotation_matrix(max(beam.get_n_max(), n_max), R=R)
        beam, D = rotate(beam, wigner=D)
        return beam, D
    elif wigner.size and not R.size:
        if False:
            pass
            #TODO Add case for multiple beams
#            ibeam = beam;
#            beam = ott.Bsc();#
#
#            for ii = 1:numel(p.Results.wigner):
#                sz = size(ibeam.a, 1);
#                D2 = p.Results.wigner{ii}(1:sz, 1:sz);
#                beam = beam.append(D2 * ibeam);        
        else:
            sz = beam.a.shape[0]
            D2 = wigner[:sz, :sz]
            beam = D2*beam
        return beam, D2
    else:
        raise ValueError('One of wigner or R must be specified')
