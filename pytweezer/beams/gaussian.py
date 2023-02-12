import numpy as np
from .point_match import PointMatch

class Gaussian(PointMatch):
    def __init__(self, beam_function='lg',
                        mode=np.array([0, 0]),
                        n_max=None, 
                        zero_rejection_level=1e-8,
                        offset=np.array([[0],[0],[0]]),
                        polarization=np.array([1, 1j]),
                        lambda_0=1,
                        power=0,
                        verbose=False,
                        translation_method='default',
                        omega=2*np.pi,
                        index_m=1,
                        na=0,
                        angle_deg=0,
                        angle=0,
                        truncation_angle=np.pi/2,
                        angular_scaling='tantheta'
                        ):
        pass

'''
        if self.validate_beam_type(beam_function):
            self.beam_function = beam_function
        else:
            raise ValueError('Beam type inserted not valid! Allowed values are \`lg\`, \`hg\` or \`ig\`.')
        self.mode = self.mode
        self.offset = self.offset
        self.polarization = polarization

        if self.validate_translation_method(translation_method):
            self.translation_method = translation_method
        else:
            raise ValueError('Translation method inserted \
                not valid! Allowed values are \`default\` or \`new_beam\`.')
        self.compute_k_medium(index_n, lambda_0)
        self.omega = omega
        self.offset = offset
        self.truncation_angle = self.pi/2
        axisymmetry = 1
        radial = 0
        azimuthal = 0

        if self.beam_function == 'hg':
            assert mode.size == 2, 'Mode must be array with two elements'
            m, n = self.mode[0], self.mode[1]
            paraxial_order = n+m
            modeweights,initial_mode,final_mode = paraxial_transformation_matrix(paraxial_order,0,1,0);
            row = find(final_mode(:,1)==m,1)
        elif self.beam_function == 'lg':

            assert(numel(p.Results.mode) == 2, ...
            'ott:BscPmGauss:wrong_mode_length', ...
            'mode must be 2 element vector');

            radial_mode = p.Results.mode(1);
            azimuthal_mode = p.Results.mode(2);
            assert(radial_mode - floor(radial_mode) == 0 && radial_mode >= 0, ...
            'ott:BscPmGauss:invalid_radial_mode', ...
            'Radial mode index must be positive integer');
            assert(azimuthal_mode - floor(azimuthal_mode) == 0, ...
            'ott:BscPmGauss:invalid_azimuthal_mode', ...
            'Azimuthal mode index must be integer');

            paraxial_order=2*radial_mode+abs(azimuthal_mode);
            modeweights=eye(paraxial_order+1);
            row=(azimuthal_mode+paraxial_order)/2+1;
            
            i2_out= (-paraxial_order:2:paraxial_order).';
            i1_out=floor((paraxial_order-abs(i2_out))/2);
            
            initial_mode=[i1_out,i2_out];
        elif self.beam_function == 'ig':
            assert(numel(p.Results.mode) == 4, ...
            'ott:BscPmGauss:wrong_mode_length', ...
            'mode must be 4 element vector');

            paraxial_order = p.Results.mode(1);
            azimuthal_mode = p.Results.mode(2);
            parity = p.Results.mode(3);
            elipticity = p.Results.mode(4);
            
            [modeweights,initial_mode,final_mode] = ...
                paraxial_transformation_matrix(paraxial_order,0,[2,elipticity],0);
            
            [row]=find(and(final_mode(:,2)==azimuthal_mode, ...
                final_mode(:,3)==parity),1);
            
            if and(paraxial_order>1,isempty(row))
            ott.warning('external');
            error('Observe parity convensions!')
            end
        keepz=(abs(modeweights(row,:))>0);
        initial_mode=initial_mode(keepz,:);
        c=modeweights(row,keepz);

        beam_angle_specified = ~isempty(p.Results.angle) ...
            + ~isempty(p.Results.angle_deg) + ~isempty(p.Results.NA);

        if beam_angle_specified > 1
            ott.warning('external');
            error('Too many inputs.  Only specify NA, angle or angle_deg');
        elseif isempty(p.Results.angle) && isempty(p.Results.angle_deg)
            if isempty(p.Results.NA)
                NA = 1.02;
                index = 1.33;
            else
                NA = p.Results.NA;
            if ~isempty(p.Results.index_medium)
                index = p.Results.index_medium;
            else
                ott.warning('external');
                error('Need to specify index_medium with NA');
          
            beam_angle_deg = asin(NA/index)*180.0/pi;
        elseif ~isempty(p.Results.angle_deg)
            beam_angle_deg = p.Results.angle_deg;
        elseif ~isempty(p.Results.angle)
            beam_angle_deg = p.Results.angle*180/pi;
        beam.angle = beam_angle_deg * pi/180;

        xcomponent = p.Results.polarisation(1);
        ycomponent = p.Results.polarisation(2);
        offset = p.Results.offset;

        if numel(offset) == 3 && any(abs(offset(1:2))>0)
        
        % Only warn if using beams that support matrix translations
            if strcmpi(p.Results.translation_method, 'Default')
                ott.warning('external');
                ott.warning('ott:bsc_pointmatch_farfield:offsets', ...
                    ['Beam offsets with x and y components cannot be ' ...
                    'axi-symmetric, beam symmetry is now off, and the ' ...
                    'calculation will be much slower. It is highly recommended ' ...
                    'that a combination of rotations and translations are ' ...
                    'used on BSCs instead.']);
                    ott.warning('internal');
        
            axisymmetry=0;
      
        w0 = ott.utils.paraxial_beam_waist(paraxial_order);
        wscaling=1/tan(abs(beam_angle_deg/180*pi));
        if isempty(p.Results.Nmax)
            nmax = 100;
        else
            nmax = p.Results.Nmax;
        ntheta = (nmax + 1);
        nphi = 2*(nmax + 1);
        if axisymmetry
            ntheta = 2*(nmax+1);
            nphi = 3;
            if ~strcmp(beam.gtype, 'lg')
                nphi = paraxial_order+3-rem(paraxial_order,2);
        if ~isempty(p.Results.offset)
            offset_lambda = vecnorm(p.Results.offset)*beam.k_medium/(2*pi);
            ntheta = max(ntheta, 3*ceil(offset_lambda));
            nphi = max(nphi, 2*3*ceil(offset_lambda));
        [theta,phi] = ott.utils.angulargrid(ntheta,nphi);
        np = length(theta);
        central_amplitude = 1;
        rw = 2*(wscaling * w0)^2 * tan(theta).^2 ;
        dr = (wscaling * w0) * (sec(theta)).^2 ;
      
      if strcmpi(p.Results.angular_scaling, 'tantheta')
        % Nothing to do
      elseif strcmpi(p.Results.angular_scaling, 'sintheta')

        wscaling=1/sin(abs(beam_angle_deg/180*pi));

        rw = 2*(wscaling * w0)^2 * sin(theta).^2 ;
        dr = (wscaling * w0) * abs(cos(theta)) ;
        
      else
        error('Unknown angular_scaling parameter value');
      end
      
      beam.angular_scaling = p.Results.angular_scaling;

      % degree and order of all modes
      total_modes = nmax^2 + 2*nmax;
      [nn,mm] = ott.utils.combined_index((1:total_modes)');

      mode_index_vector=[];
      beam_envelope = zeros(np,length(c));
      for ii=1:length(c)
          radial_mode=initial_mode(ii,1);
          azimuthal_mode=initial_mode(ii,2);
          
          norm_paraxial=sqrt(2*factorial(radial_mode)/(pi*factorial(radial_mode+abs(azimuthal_mode))));
          L = laguerre(radial_mode,abs(azimuthal_mode),rw);
          beam_envelope(:,ii) = norm_paraxial.*rw.^abs(azimuthal_mode/2) .* L .* exp(-rw/2 + 1i*azimuthal_mode*phi+1i*pi/2*(radial_mode*2+abs(azimuthal_mode)+1));
          mode_input_power=sqrt(sum(2*pi*abs(beam_envelope(:,ii)).^2.*sqrt(rw/2).*abs(dr)));
          aperture_power_normalization=sqrt(sum(2*pi*abs(beam_envelope(:,ii)).^2.*sin(theta)));
          
          beam_envelope(:,ii)=c(ii)*beam_envelope(:,ii)/aperture_power_normalization*mode_input_power;
          
          mode_index_vector=[mode_index_vector; ...
              find(mm==azimuthal_mode+1-max([azimuthal,radial]) ...
              | mm==azimuthal_mode-1+max([azimuthal,radial]))];

      end
      mode_index_vector=unique(mode_index_vector);

      beam_envelope=sum(beam_envelope,2);
      outbeam = theta < pi-beam.truncation_angle;
      beam_envelope(outbeam) = 0;

      if ~isempty(offset)
        rhat = rtpv2xyzv( ones(size(theta)), zeros(size(theta)), ...
            zeros(size(theta)), ones(size(theta)), theta, phi );
        [offset,rhat] = matchsize(offset.',rhat);
        phase_shift = exp( 1i * beam.k_medium * dot(offset,rhat,2) );
        beam_envelope = beam_envelope .* phase_shift;
      end
      Ex = xcomponent * beam_envelope * central_amplitude;
      Ey = ycomponent * beam_envelope * central_amplitude;

      if any(azimuthal|radial)
        Etheta=-radial*xcomponent*beam_envelope * central_amplitude;
        Ephi=azimuthal*ycomponent*beam_envelope * central_amplitude;
      else
        Etheta = - Ex .* cos(phi) - Ey .* sin(phi);
        Ephi = - Ex .* sin(phi) + Ey .* cos(phi);
      end

      e_field = [ Etheta(:); Ephi(:) ];

      if axisymmetry
        nn=nn(mode_index_vector);
        mm=mm(mode_index_vector);

        removeels=find(abs(mm)>paraxial_order+1);
        nn(removeels)=[];
        mm(removeels)=[];
      end

      % Do the point matching and store the result
      [beam.a, beam.b] = beam.bsc_farfield(nn, mm, e_field, theta, phi, ...
        'zero_rejection_level', p.Results.zero_rejection_level);

      % If no Nmax supplied, shrink the beam to the smallest size that
      % preserves the beam power
      if isempty(p.Results.Nmax)
        beam = beam.shrink_Nmax();
      end

      % Normalize the beam power
      if ~isempty(p.Results.power)
        beam.power = p.Results.power;
      end

      ott.warning('external');
    end

      #  self.beam = PointMatch(**kwargs)
        self.beam_type = 'incident'
        self.beam.basis = 'regular'

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

    @staticmethod
    def supported_beam_type(self, s):
        return True if s in ('lg', 'hg', 'ig') else False

    @staticmethod
    def validate_translation_method(self, method):
        return True if method in ('default', 'new_beam') else False


'''

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