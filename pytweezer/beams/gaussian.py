import numpy as np
from .point_match import PointMatch
from pytweezer.utils import angular_grid, combined_index, laguerre, paraxial_beam_waist
from numpy.linalg import norm

class Gaussian(PointMatch):
    def __init__(self, beam_function='lg',
                        mode=np.array([0, 0]),
                        n_max=100, 
                        zero_rejection_level=1e-8,
                        offset=np.array([[0],[0],[0]]),
                        polarization=np.array([1, 1j]),
                        lambda_0=1,
                        power=0,
                        verbose=False,
                        translation_method='default',
                        omega=2*np.pi,
                        index_m=1.33,
                        na=1.02,
                        angle=0,
                        truncation_angle=np.pi/2,
                        angular_scaling='tantheta'
                        ):
        if self.validate_beam_type(beam_function):
            self.beam_function = beam_function
        else:
            raise ValueError('Beam type inserted not valid! Allowed values are \`lg\`, \`hg\` or \`ig\`.')
        self.mode = mode
        self.offset = offset
        self.polarization = polarization

        if self.validate_translation_method(translation_method):
            self.translation_method = translation_method
        else:
            raise ValueError('Translation method inserted \
                not valid! Allowed values are \`default\` or \`new_beam\`.')
        self.compute_k_medium(index_m, lambda_0)
        self.omega = omega
        self.offset = offset
        self.truncation_angle = np.pi/2
        axi_symmetry, radial, azimuthal = 1, 0, 0
        if self.beam_function == 'hg':
            raise ValueError('Beam function not implemented yet!')
        elif self.beam_function == 'lg':
            assert mode.size == 2, 'Wrong mode length, for LG mode must be array containig 2 elements'
            radial_mode = self.mode[0]
            azimuthal_mode = self.mode[1]
            assert (radial_mode - np.floor(radial_mode)) == 0 and radial_mode >= 0, 'Radial mode must be positive integer, condition not attended'
            assert (azimuthal_mode - np.floor(azimuthal_mode)) == 0, 'Azimuthal mode index must be integer'
            paraxial_order = 2*radial_mode + np.abs(azimuthal_mode)
            mode_weights = np.eye(paraxial_order+1)
            row = (azimuthal_mode+paraxial_order)/2+1            
            i2_out = np.arange(-paraxial_order, paraxial_order+1, 2).T
            i1_out = np.floor((paraxial_order-np.abs(i2_out))/2)
            initial_mode = np.array([i1_out, i2_out]).T
        elif self.beam_function == 'ig':
            raise ValueError('Beam function not implemented yet!')
        keepz = np.argwhere(np.abs(mode_weights[int(row-1),:]) > 0)
        initial_mode = initial_mode[keepz, :][0]
        c = mode_weights[int(row-1), keepz]
        self.angle = np.arcsin(na/self.index_m)
        x_comp = polarization[0]
        y_comp = polarization[1]
        if offset.size == 3 and (np.abs(offset[:2])>0).any():
            if self.translation_method == 'default':
                warnings.warn('Beam offsets with x and y components cannot be \
                    axi-symmetric, beam symmetry is now off, and the \
                    calculation will be much slower. It is highly recommended \
                    that a combination of rotations and translations are \
                    used on BSCs instead.')                
            axi_symmetry = 0
        w0 = paraxial_beam_waist(paraxial_order)
        wscaling = 1/np.tan(np.abs(self.angle))
        n_theta = (n_max + 1)
        n_phi = 2*(n_max + 1)
        if axi_symmetry:
            n_theta = 2*(n_max+1)
            n_phi = 3
            if self.beam_function == 'lg':
                n_phi = paraxial_order + 3 - np.remainder(paraxial_order,2)
        if norm(self.offset - np.array([[0], [0], [0]])) > 1e-12:
            offset_lambda = norm(offset)*self.k_medium/(2*np.pi)
            n_theta = max(n_theta, 3*np.ceil(offset_lambda))
            n_phi = max(n_phi, 2*3*np.ceil(offset_lambda))
        theta, phi = angular_grid(n_theta, n_phi)
        n_p = theta.size
        central_amplitude = 1
        rw = 2*np.power(wscaling * w0, 2) * np.power(np.tan(theta),2)
        dr = (wscaling * w0) * np.power(1/np.cos(theta),2)
        if angular_scaling == 'tantheta':
            pass
        elif angular_scaling == 'sintheta':
            wscaling = 1/np.sin(np.abs(self.angle))
            rw = 2*np.power(wscaling * w0, 2) * np.power(np.sin(theta), 2)
            dr = (wscaling * w0)*np.abs(np.cos(theta))
        else:
            raise ValueError('Unknown angular_scaling parameter value')
        self.angular_scaling = angular_scaling
        total_modes = np.power(n_max, 2) + 2*n_max
        nn, mm = combined_index(np.arange(1, total_modes+1, 1).T)
        nn, mm = nn.astype(int), mm.astype(int)
        mode_index_vector = []
        beam_envelope = np.zeros((n_p, c.size), dtype=complex)
        for i in range(1, c.size+1):
            radial_mode = initial_mode[i-1, 0]
            azimuthal_mode = initial_mode[i-1,1]            
            norm_paraxial = np.sqrt(2*np.math.factorial(radial_mode)/(np.pi*np.math.factorial(radial_mode+np.abs(azimuthal_mode))))
            L = laguerre(int(radial_mode), int(np.abs(azimuthal_mode)), rw)
            mult1 = np.exp(-rw/2 + 1j*azimuthal_mode*phi+1j*np.pi/2*(radial_mode*2+np.abs(azimuthal_mode)+1))
            elements = norm_paraxial*np.power(rw, np.abs(azimuthal_mode/2))*L*mult1
            
            beam_envelope[:,i-1] = norm_paraxial*np.power(rw, np.abs(azimuthal_mode/2))*L*mult1
            mode_input_power = np.sqrt((2*np.pi*np.power(abs(beam_envelope[:, i-1]),2)*np.sqrt(rw/2)*np.abs(dr)).sum())
            aperture_power_normalization = np.sqrt(sum(2*np.pi*np.power(np.abs(beam_envelope[:,i-1]), 2)*np.sin(theta)))
            beam_envelope[:,i-1] = c[i-1]*beam_envelope[:, i-1]/aperture_power_normalization*mode_input_power
            max_aux = max([azimuthal,radial]) 
            conditional_args = np.where((mm==azimuthal_mode+1-max_aux) | (mm==azimuthal_mode-1+max_aux))
            if conditional_args[0].size:
                conditional_args = tuple(conditional_args[0].T)
                mode_index_vector = list(conditional_args)+mode_index_vector
        mode_index_vector = set(mode_index_vector)
        beam_envelope = beam_envelope.sum(axis=1)
        outbeam = theta < np.pi - self.truncation_angle
        beam_envelope[outbeam] = 0
        if not np.array_equal(offset, np.array([[0], [0], [0]])):
            rhat = rtpv2xyzv( np.ones(size(theta)), zeros(size(theta)),
                np.zeros(size(theta)), np.ones(size(theta)), theta, phi)
            offset, rhat = match_size(offset.T,rhat)
            phase_shift = np.exp( 1j * beam.k_medium * dot(offset,rhat,2))
            beam_envelope = beam_envelope * phase_shift
        Ex = x_comp * beam_envelope * central_amplitude
        Ey = y_comp * beam_envelope * central_amplitude

        if azimuthal or radial:
            E_theta = -radial*x_comp * beam_envelope * central_amplitude;
            E_phi = azimuthal*y_comp * beam_envelope * central_amplitude;
        else:
            E_theta = - Ex * np.cos(phi) - Ey * np.sin(phi)
            E_phi = - Ex * np.sin(phi) + Ey * np.cos(phi)
        e_field = np.concatenate([E_theta[:], E_phi[:]])
        if axi_symmetry:
            nn = nn[list(mode_index_vector)]
            mm = mm[list(mode_index_vector)]
            mm = mm[abs(mm) <= paraxial_order+1]
            nn = nn[abs(mm) <= paraxial_order+1]
            nn = np.sort(nn)
        a, b, _ = self.bsc_far_field(nn, mm, e_field, theta, phi, zero_rejection_level= zero_rejection_level)
        self.a = a
        self.b = b
        self.power = power
        self.beam_type = 'incident'
        self.beam_basis = 'regular'

    def translate_z(self, **kwargs):
        if self.translation_method == 'default':        
            A, B = self.translate_z_point_match(**kwargs)
            return _, A, B
        elif self.translation_method == 'new_beam_offset':
            if not kwargs.get('z'):
                kwargs['z'] = 0
            if not kwargs.get('n_max'):
                kwargs['n_max'] = self.n_max
            return Gaussian(self.gtype, self.mode, 
                offset=self.offset+np.array([0,0,kwargs['z']]),
                omega=self.omega, power=self.power, lambda_m=self.lambda_b,
                polarization=self.polarization, truncation_angle=self.truncation_angle,
                n_max=self.n_max, angle=self.angle)

    def validate_beam_type(self, s):
        return True if s in ('lg', 'hg', 'ig') else False

    def validate_translation_method(self, method):
        return True if method in ('default', 'new_beam') else False


'''


    
    function varargout = translateXyz(beam, varargin)
      %TRANSLATEXYZ translate the beam given Cartesian coordinates
      %
      % beam = TRANSLATEXYZ(xyz) translate the beam to locations given by
      % the xyz coordinates, where xyz is a 3xN matrix of coordinates.
      % If the translation_method is NewBeamOffset, a new beam is generated.
      %
      % TRANSLATEXYZ(Az, Bz, D) translate the beam using
      % z-translation and rotation matricies.
      %
      % [beam, Az, Bz, D] = TRANSLATEXYZ(...) returns the z-translation
      % matrices, the rotation matrix D, and the translated beam.
      %
      % [beam, A, B] = TRANSLATEXYZ(...) returns the translation matrices
      % and the translated beam.
      %
      % [beam, AB] = TRANSLATEXYZ(...) returns the A, B matricies packed
      % so they can be directly applied to the beam: tbeam = AB * beam.
      %
      % TRANSLATEXYZ(..., 'Nmax', Nmax) specifies the output beam Nmax.
      % Takes advantage of not needing to calculate a full translation matrix.

      if strcmpi(beam.translation_method, 'Default')
        
        % Use translation matrix method
        [varargout{1:nargout}] = translateXyz@ott.BscPointmatch(beam, varargin{:});
        
      elseif strcmpi(beam.translation_method, 'NewBeamOffset')
        
        p = inputParser;
        p.addOptional('opt1', []);    % xyz or Az
        p.addOptional('opt2', []);    % [] or Bz
        p.addOptional('opt3', []);    % [] or D
        p.addParameter('Nmax', beam.Nmax);
        p.parse(varargin{:});
        
        assert(isempty(p.Results.opt2) && isempty(p.Results.opt3), ...
          'Rotation and translation matries not supported with this method');
        
        % Generate the new beam
        varargout{1} = ott.BscPmGauss(beam.gtype, beam.mode, ...
          'offset', beam.offset + p.Results.opt1, ...
          'omega', beam.omega, 'power', beam.power, ...
          'wavelength_medium', beam.wavelength, ...
          'polarisation', beam.polarisation, ...
          'truncation_angle', beam.truncation_angle, ...
          'Nmax', p.Results.Nmax, 'angle', beam.angle);
        varargout{1}.type = beam.type;
        varargout{1}.basis = beam.basis;
      end
    end
    
    function varargout = translateRtp(beam, varargin)
      %TRANSLATERTP translate the beam given spherical coordinates
      %
      % beam = TRANSLATERTP(rtp) translate the beam to locations given by
      % the xyz coordinates, where rtp is a 3xN matrix of coordinates.
      % If the translation_method is NewBeamOffset, a new beam is generated.
      %
      % TRANSLATERTP(Az, Bz, D) translate the beam using
      % z-translation and rotation matricies.
      %
      % [beam, Az, Bz, D] = TRANSLATERTP(...) returns the z-translation
      % matrices, the rotation matrix D, and the translated beam.
      %
      % [beam, A, B] = TRANSLATERTP(...) returns the translation matrices
      % and the translated beam.
      %
      % [beam, AB] = TRANSLATERTP(...) returns the A, B matricies packed
      % so they can be directly applied to the beam: tbeam = AB * beam.
      %
      % TRANSLATERTP(..., 'Nmax', Nmax) specifies the output beam Nmax.
      % Takes advantage of not needing to calculate a full translation matrix.
      
      if strcmpi(beam.translation_method, 'Default')
        
        % Use translation matrix method
        [varargout{1:nargout}] = translateRtp@ott.BscPointmatch(beam, varargin{:});
        
      elseif strcmpi(beam.translation_method, 'NewBeamOffset')
        
        p = inputParser;
        p.addOptional('opt1', []);    % xyz or Az
        p.addOptional('opt2', []);    % [] or Bz
        p.addOptional('opt3', []);    % [] or D
        p.addParameter('Nmax', beam.Nmax);
        p.parse(varargin{:});
        
        assert(isempty(p.Results.opt2) && isempty(p.Results.opt3), ...
          'Rotation and translation matries not supported with this method');
        
        % Generate the new beam
        xyz = ott.utils.rtp2xyz(p.Results.opt1);
        varargout{1} = ott.BscPmGauss(beam.gtype, beam.mode, ...
          'offset', beam.offset + xyz(:), ...
          'omega', beam.omega, 'power', beam.power, ...
          'wavelength_medium', beam.wavelength, ...
          'polarisation', beam.polarisation, ...
          'truncation_angle', beam.truncation_angle, ...
          'Nmax', p.Results.Nmax, 'angle', beam.angle);
        varargout{1}.type = beam.type;
        varargout{1}.basis = beam.basis;
      end
    end
  end
end

'''