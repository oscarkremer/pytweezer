def translate(beam, A, B):
    beam = [ A, B, B, A ] * beam;
    return beam


def _translate_z(beam, z=0):
beam.translation_method == 'default'"

        if not kwargs.get('z'):
            raise ValueError('No z informed for translation method!')
        else:
            z = kwargs['z']
            beam.dz = beam.dz + np.abs(z)
            z = z * beam.k_m / 2 / np.pi
            i_beam = beam
            if isinstance(z, np.ndarray):
                for i in range(1, z.size+1):
                    A, B = ibeam.translateZ_type_helper(z[i], [p.Results.Nmax, ibeam.Nmax]);
                    beam = beam.append(ibeam.translate(A, B));
                    beam.basis = 'regular';
            else:
                A, B = ibeam.translateZ_type_helper(z, [p.Results.Nmax, ibeam.Nmax]);
                beam = beam.append(ibeam.translate(A, B))
                beam.basis = 'regular';
        return beam, A, B


def translate_z(beam, **kwargs):
    if beam.translation_method == 'default':        
        A, B = _translate_z(beam, **kwargs)
        return _, A, B
    elif beam.translation_method == 'new_beam_offset':
        if not kwargs.get('z'):
            kwargs['z'] = 0
        if not kwargs.get('n_max'):
            kwargs['n_max'] = self.n_max
        return Gaussian(beam.gtype, self.mode, 
            offset=self.offset+np.array([0,0,kwargs['z']]),
            omega=self.omega, power=self.power, lambda_m=self.lambda_b,
            polarization=self.polarization, truncation_angle=self.truncation_angle,
            n_max=self.n_max, angle=self.angle)

