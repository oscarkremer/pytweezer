import numpy as np
from copy import copy
from pytweezer.utils import translate_z as translation, combined_index, xyz2rtp
from .rotation import rotate_yz, rotate
import time

def translate(beam, A, B):
    AB = np.block([[A, B], [B, A]])
    return copy(beam)*AB

def _translate_z_(beam, z=0, n_max=None):
    if not n_max:
        n_max = beam.get_n_max()
    def translate_z_type_helper(z, n_max, basis):
        hash_basis_type = {'incoming': 'sbesselh2', 'outgoing': 'sbesselh1', 'regular': 'sbesselj'}
        translation_type = hash_basis_type[basis]
        A, B, _ = translation(n_max, z, function_type=translation_type)
        return A, B
    if isinstance(z, np.ndarray):
        z = z[0] if z.size == 1 else z
    beam.dz = beam.dz + np.abs(z)
    z = z * beam.k_m / 2 / np.pi
    if isinstance(z, np.ndarray):
        for i in range(z.size):
            A, B = translate_z_type_helper(z[i], 
                np.array([n_max, beam.get_n_max()]).astype(int), 
                beam.get_basis())
            beam.append(translate(beam, A, B))
            beam.set_basis('regular')
    else:
        A, B = translate_z_type_helper(z, 
            np.array([n_max, beam.get_n_max()]).astype(int), 
            beam.get_basis())
        beam.append(translate(beam, A, B))
        beam.set_basis('regular')
    return beam, A, B

def translate_xyz(beam, position, n_max=100):
    if len(position.shape)==2:
        position = position.T
        rtp = xyz2rtp(x=position[:,0],y=position[:,1], z=position[:,2]).T
    elif len(position.shape) == 1:
        r, t, p = xyz2rtp(x=position[0],y=position[1], z=position[2])
        rtp = np.array([[r, t, p]]).T
    end = time.time()
    return translate_rtp(beam, rtp, n_max=n_max)[0]
    
def translate_rtp(beam, position, n_max=100):
    if n_max.size == 1:
        o_n_max = n_max
    elif n_max.size == 2:
        o_n_max = n_max[1]
    else:
        raise ValueError('Nmax must be 2 element vector or scalar');
    r = position[0, :]
    theta = position[1, :]
    phi = position[2, :]
    if any((theta != 0 and abs(theta) != np.pi) or phi != 0):
        i_beam = copy(beam)
        for i in range(r.size):
            tbeam, D = rotate_yz(i_beam, theta[i], phi[i], r[i], n_max= max(o_n_max, beam.n_max))
            tbeam, A, B = translate_z(tbeam, r[i], n_max=o_n_max)
            tbeam, _ = rotate(tbeam, wigner=D.T)
            beam.append(tbeam)
        return beam, A, B, D
    else:
        d_n_max = max(o_n_max, beam.n_max)
#        print(combined_index(d_n_max, d_n_max), d_n_max)
        D = np.eye(combined_index(d_n_max, d_n_max))
        idx = abs(theta) == np.pi
#        print(idx, np.where(theta==))
        r[idx] = -r[idx]
        start = time.time()
        if r.size == 1:
            beam, A, B = translate_z(beam, r, n_max=o_n_max)
        else:
            beam, A, B = translate_z(beam, r, n_max=o_n_max)
        end = time.time()
        return beam, A, B, D


def translate_z(beam, z, n_max=None):
    if beam.translation_method == 'default':        
        beam, A, B = _translate_z_(beam, z=z, n_max=n_max)
        return beam, A, B
    elif beam.translation_method == 'new_beam_offset':
        if not kwargs.get('z'):
            kwargs['z'] = 0
        if not kwargs.get('n_max'):
            n_max = beam.n_max
        return Gaussian(beam.gtype, beam.mode, 
            offset=beam.offset+np.array([0,0,z]),
            omega=beam.omega, power=beam.power, lambda_m=beam.lambda_b,
            polarization=beam.polarization, truncation_angle=beam.truncation_angle,
            n_max=n_max, angle=beam.angle), A, B

