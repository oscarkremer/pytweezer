import numpy as np
from copy import copy
from pytweezer.utils import translate_z as translation

def translate(beam, A, B):
    AB = np.block([[A, B], [B, A]])
    return copy(beam)*AB

def _translate_z_(beam, z=0):
    def translate_z_type_helper(z, n_max, basis):
        if basis == 'incoming':
            translation_type = 'sbesselh2';
        elif basis == 'outgoing':
            translation_type = 'sbesselh1';
        elif basis == 'regular':
            translation_type = 'sbesselj'
        A, B, _ = translation(n_max, z, function_type=translation_type)
        return A, B
    beam.dz = beam.dz + np.abs(z)
    z = z * beam.k_m / 2 / np.pi
    print(z)
    if isinstance(z, np.ndarray):
        for i in range(1, z.size+1):
            A, B = translate_z_type_helper(z[i], beam.n_max)
            print(A.sum())
            beam.append(translate(beam, A, B))
            beam.beam_basis = 'regular'
    else:
        print(beam.n_max)
        A, B = translate_z_type_helper(z, beam.n_max, beam.beam_basis)
        beam.append(translate(beam, A, B))
        beam.beam_basis = 'regular'
    return beam, A, B


def translate_z(beam, z, n_max=None):
    if beam.translation_method == 'default':        
        beam, A, B = _translate_z_(beam, z=z)
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

