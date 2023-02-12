import warnings
import numpy as np
from pytweezer.utils import combined_index, translate_z

class Beam:
    def __init__(self, a, b, basis, beam_type, k_m=2*np.pi, omega=2*np.pi, dz=0):
        self.dz = dz
        self.k_medium = k_m
        self.omega = omega

        assert a.shape==b.shape, 'Mismatch between the shapes of the coefficients'      
        if isinstance(a, np.ndarray):
            if len(a.shape) = 1
                a = a.reshape((a.shape[0], 1))
        if isinstance(b, np.ndarray):
            if len(b.shape) = 1
                b = b.reshape((b.shape[0], 1))
        assert a.shape[0] >= 3, 'Number of multiploe must be at least 3'
        assert np.sqrt(a.shape[0]+1) == np.floor(np.sqrt(a.shape[0]+1)), 
            'Number of multipoles must be 3, 8, 15, 24, ...'
        self.a = a
        self.b = b
        self.basis = basis
        self.type = beam_type


'''
 function data = GetVisualisationData(field_type, xyz, rtp, vxyz, vrtp)
      % Helper to generate the visualisation data output.
      % This function is not intended to be called directly, instead
      % see :meth:`visualise` or :meth:`visualiseFarfield`.
      %
      % Usage
      %   GetVisualisationData(field_type, xyz, rtp, vxyz, vrtp)
      %   Takes a field_type string, the coordinates (either xyz or rtp),
      %   and the data values (either vxyz or vrtp).
      %
      % Parameters
      %   - xyz, rtp, vxyz, vrtp -- (Nx3 numeric) Coordinates in a
      %     suitable form to be passed to `ott.utils.xyz2rtp` and similar
      %     functions.  Pass empty arrays for unused values.
      %   - field_type -- (enum) Type of field to calculate.
      %     Supported types include:
      %       - 'irradiance'  -- :math:`\sqrt{|Ex|^2 + |Ey|^2 + |Ez|^2}`
      %       - 'E2' -- :math:`|Ex|^2 + |Ey|^2 + |Ez|^2`
      %       - 'Sum(Abs(E))' -- :math:`|Ex| + |Ey| + |Ez|`
      %
      %       - Re(Er), Re(Et), Re(Ep), Re(Ex), Re(Ey), Re(Ez)
      %       - Abs(Er), Abs(Et), Abs(Ep), Abs(Ex), Abs(Ey), Abs(Ez)
      %       - Arg(Er), Arg(Et), Arg(Ep), Arg(Ex), Arg(Ey), Arg(Ez)
      
      assert(size(xyz, 2) == 3 || size(xyz, 2) == 0, ...
        'xyz must be Nx3 matrix');
      assert(size(vxyz, 2) == 3 || size(vxyz, 2) == 0, ...
        'vxyz must be Nx3 matrix');
      assert(size(rtp, 2) == 3 || size(rtp, 2) == 0, ...
        'rtp must be Nx3 matrix');
      assert(size(vrtp, 2) == 3 || size(vrtp, 2) == 0, ...
        'vrtp must be Nx3 matrix');
      
      % Get the coordinates
      if isempty(xyz) && ~isempty(rtp)
        xyz = ott.utils.rtp2xyz(rtp);
      elseif isempty(rtp) && ~isempty(xyz)
        rtp = ott.utils.xyz2rtp(xyz);
      elseif isempty(rpt) && isempty(xyz)
        error('OTT:BSC:GetVisualisationData:no_coords', ...
          'Must supply coordinates');
      end
      
      % Get the data
      if isempty(vxyz) && ~isempty(vrtp)
        vxyz = ott.utils.rtpv2xyzv(vrtp, rtp);
      elseif isempty(vrtp) && ~isempty(vxyz)
        vrtp = ott.utils.xyzv2rtpv(vxyz, xyz);
      elseif isempty(vrtp) && isempty(vxyz)
        error('OTT:BSC:GetVisualisationData:no_data', ...
          'Must supply data');
      else
        error('OTT:BSC:GetVisualisationData:too_much_data', ...
          'Must supply only one data variable');
      end
      
      % Generate the requested field
      if strcmpi(field_type, 'irradiance')
        data = sqrt(sum(abs(vxyz).^2, 2));
      elseif strcmpi(field_type, 'E2')
        data = sum(abs(vxyz).^2, 2);
      elseif strcmpi(field_type, 'Sum(Abs(E))')
        data = sum(abs(vxyz), 2);
        
      elseif strcmpi(field_type, 'Re(Er)')
        data = real(vrtp(:, 1));
      elseif strcmpi(field_type, 'Re(Et)')
        data = real(vrtp(:, 2));
      elseif strcmpi(field_type, 'Re(Ep)')
        data = real(vrtp(:, 3));
        
      elseif strcmpi(field_type, 'Re(Ex)')
        data = real(vxyz(:, 1));
      elseif strcmpi(field_type, 'Re(Ey)')
        data = real(vxyz(:, 2));
      elseif strcmpi(field_type, 'Re(Ez)')
        data = real(vxyz(:, 3));
        
      elseif strcmpi(field_type, 'Abs(Er)')
        data = abs(vrtp(:, 1));
      elseif strcmpi(field_type, 'Abs(Et)')
        data = abs(vrtp(:, 2));
      elseif strcmpi(field_type, 'Abs(Ep)')
        data = abs(vrtp(:, 3));
        
      elseif strcmpi(field_type, 'Abs(Ex)')
        data = abs(vxyz(:, 1));
      elseif strcmpi(field_type, 'Abs(Ey)')
        data = abs(vxyz(:, 2));
      elseif strcmpi(field_type, 'Abs(Ez)')
        data = abs(vxyz(:, 3));
        
      elseif strcmpi(field_type, 'Arg(Er)')
        data = angle(vrtp(:, 1));
      elseif strcmpi(field_type, 'Arg(Et)')
        data = angle(vrtp(:, 2));
      elseif strcmpi(field_type, 'Arg(Ep)')
        data = angle(vrtp(:, 3));
        
      elseif strcmpi(field_type, 'Arg(Ex)')
        data = angle(vxyz(:, 1));
      elseif strcmpi(field_type, 'Arg(Ey)')
        data = angle(vxyz(:, 2));
      elseif strcmpi(field_type, 'Arg(Ez)')
        data = angle(vxyz(:, 3));
        
      elseif strcmpi(field_type, 'Er')
        data = vrtp(:, 1);
      elseif strcmpi(field_type, 'Et')
        data = vrtp(:, 2);
      elseif strcmpi(field_type, 'Ep')
        data = vrtp(:, 3);
        
      elseif strcmpi(field_type, 'Ex')
        data = vxyz(:, 1);
      elseif strcmpi(field_type, 'Ey')
        data = vxyz(:, 2);
      elseif strcmpi(field_type, 'Ez')
        data = vxyz(:, 3);

      else
        error('OTT:BSC:GetVisualisationData:unknown_field_type', ...
          'Unknown field type value');
      end
      
    end
  end
'''
    @staticmethod
    def make_beam_vector(self, a, b, n, m, n_max):
        if not n.shape or not n:
            warnings.warn('No modes in beam')
        if not n.shape or not n:
            n_max = 0
        else:
            n_max = n.max()

        total_orders = combined_index(n_max, n_max)
        ci = combined_index(n, m)
        nbeams = a.shape[1]
        ci, c_in_beams = np.meshgrid(ci, np.arange(1, n_beams+1, 1), indexing='ij')
        a = csr_matrix((ci, c_in_beams, a), shape=(total_orders, n_beams)).toarray()
        a = csr_matrix((ci, c_in_beams, b), shape=(total_orders, n_beams)).toarray()        
        n, m = combined_index(np.arange(1, n_max**2+2*n_max+1, 1))
        n = n.T
        m = m.T
        return a, b, n, m


    @staticmethod
    def parser_k_medium(self, **kwargs):
        if kwargs.get('k_m') and kwargs.get('lambda_m'):
            self.k_m = kwargs.get('k_m')
            warnings.warn('Both k_m and lambda_m defined, only k_m is going to \
                be used.')
        elif kwargs.get('k_m'):
            self.k_m = kwargs.get('k_m')
        elif kwargs.get('lambda_m'):
            self.k_m = 2.0*np.pi/kwargs['lambda_m']
        elif kwargs.get('index_m'):
            if kwargs.get('lambda_0'):   
                self.k_m = kwargs.get('index_m')*2*np.pi/kwargs.get('lambda_0')
            else:
                raise ValueError('Wavelength for vacuum must be specified to use medium refraction index')
        else:
            raise ValueError('Unable to determine k_medium from inputs')
    
    def translate_z_type_helper(self, z, n_max):
        if self.basis == 'incoming':
            translation_type = 'sbesselh2';
        elif self.basis == 'outgoing':
            translation_type = 'sbesselh1';
        elif self.basis == 'regular':
            translation_type = 'sbesselj'
        A, B = translate_z(n_max, z, function_type=translation_type)
        return A, B

    def append(self, other):
        if self.n_beams == 0:
            # DANGER: Have to implement this section
            beam = other
        else
            self.n_max = max(self.n_max, other.n_max)
            other.n_max = self.n_max
            self.a = [self.a, other.a]
            self.b = [self.b, other.b]


'''
    def paraxial_far_field(self, varargin)
        p = inputParser;
        p.addParameter('calcE', true);
        p.addParameter('calcH', nargout >= 2);
        p.addParameter('saveData', false);
        p.addParameter('data', []);
        p.addParameter('mapping', 'sin');
        p.addParameter('size', [50, 50]);
        p.addParameter('thetaMax', pi/2);
        p.addParameter('direction', 'pos');
        p.parse(varargin{:});
      
        xrange = linspace(-1, 1, p.Results.size(1));
        yrange = linspace(-1, 1, p.Results.size(2));
        [xx, yy] = meshgrid(xrange, yrange);
        phi = atan2(yy, xx);
        rr = sqrt(xx.^2 + yy.^2);
        switch p.Results.mapping
        case 'sin'
            theta = asin(rr);
        case 'tan'
            theta = atan(rr);
        case 'theta'
            theta = rr;
        otherwise
            error('Unknown mapping argument value, must be sin, tan or theta');
        end
        thetaMax = Inf;
        if ~isempty(p.Results.thetaMax)
            thetaMax = p.Results.thetaMax;
        end
        pinside = imag(theta) == 0 & theta < thetaMax;
        iphi = phi(pinside);
        itheta = theta(pinside);
      
        if strcmpi(p.Results.direction, 'neg')
            itheta = pi - itheta;
        elseif ~strcmpi(p.Results.direction, 'pos')
            error('Direction must be \'pos\' or \'neg\'');
        end
        [E, H, data] = beam.farfield(itheta(:), iphi(:), ...
            'saveData', p.Results.saveData, 'data', p.Results.data, ...
            'calcE', p.Results.calcE, 'calcH', p.Results.calcH);
        
        if nargout >= 1
        
        if p.Results.calcE
          % Generate the requested field
          dEt = beam.GetVisualisationData('Et', [], ...
            [itheta, iphi, ones(size(iphi))], [], E.');
          dEp = beam.GetVisualisationData('Ep', [], ...
            [itheta, iphi, ones(size(iphi))], [], E.');

          Et = zeros(p.Results.size);
          Et(pinside) = dEt;
          Ep = zeros(p.Results.size);
          Ep(pinside) = dEp;

          varargout{1} = Et;
          varargout{1}(:, :, 2) = Ep;
        end
        
        if nargout >= 2 && p.Results.calcH
          % Generate the requested field
          dHt = beam.GetVisualisationData('Et', [], ...
            [itheta, iphi, ones(size(iphi))], [], H.');
          dHp = beam.GetVisualisationData('Ep', [], ...
            [itheta, iphi, ones(size(iphi))], [], H.');

          Ht = zeros(p.Results.size);
          Ht(pinside) = dHt;
          Hp = zeros(p.Results.size);
          Hp(pinside) = dHp;

          varargout{2} = Ht;
          varargout{2}(:, :, 2) = Hp;
            return  varargout
        end
      end
      
      if nargout >= 3
        varargout{3} = data;
      end
    end

    def [E, H, data] = farfield(beam, theta, phi, varargin)
      ip = inputParser;
      ip.addParameter('calcE', true);
      ip.addParameter('calcH', nargout >= 2);
      ip.addParameter('saveData', false);
      ip.addParameter('data', []);
      ip.parse(varargin{:});

      [theta,phi] = ott.utils.matchsize(theta,phi);

      [theta_new,~,indY]=unique(theta);
      [phi_new,~,indP]=unique(phi);

      Etheta=zeros(length(theta),1);
      Ephi=zeros(length(theta),1);

      Htheta=zeros(length(theta),1);
      Hphi=zeros(length(theta),1);

      if strcmp(beam.basis, 'incoming')

        a = beam.a;
        b = beam.b;
        p = zeros(size(beam.a));
        q = zeros(size(beam.b));

      elseif strcmp(beam.basis, 'outgoing')

        a = zeros(size(beam.a));
        b = zeros(size(beam.a));
        p = beam.a;
        q = beam.b;

      else

        error('Regular wavefunctions go to zero in far-field');

      end

      a = ott.utils.threewide(a);
      b = ott.utils.threewide(b);
      p = ott.utils.threewide(p);
      q = ott.utils.threewide(q);

      [n,m]=ott.utils.combined_index(find(abs(beam.a)|abs(beam.b)));

      % Alocate memory for output data
      data = [];
      if ip.Results.saveData
        data = zeros(numel(indY), 0);
      end
      
      % Start a counter for accessing the data
      if ~isempty(ip.Results.data)
        dataCount = 0;
      end
      
      for nn = 1:max(n)

        vv=find(n==nn);
        if isempty(vv)
          continue;
        end

        %this makes the vectors go down in m for n.
        % has no effect if old version code.
        Nn = 1/sqrt(nn*(nn+1));

        % Create index arrays for a, b, q, p
        index=nn*(nn+1)+m(vv);
        aidx = full(a(index));
        bidx = full(b(index));
        pidx = full(p(index));
        qidx = full(q(index));
        
        if isempty(ip.Results.data)

          [~,Ytheta,Yphi] = ott.utils.spharm(nn,m(vv), ...
              theta_new,zeros(size(theta_new)));

          [PHI,M]=ndgrid(phi_new, m(vv));

          expimphi=exp(1i*M.*PHI);

          % Create full matrices (opt, R2018a)
          YthetaExpf = Ytheta(indY, :).*expimphi(indP, :);
          YphiExpf = Yphi(indY, :).*expimphi(indP, :);
          
          % Save the data if requested
          if ip.Results.saveData
            data(:, end+(1:size(Ytheta, 2))) = YthetaExpf;
            data(:, end+(1:size(Ytheta, 2))) = YphiExpf;
          end
          
        else
          
          % Load the data if present
          YthetaExpf = ip.Results.data(:, dataCount+(1:length(vv)));
          dataCount = dataCount + length(vv);
          YphiExpf = ip.Results.data(:, dataCount+(1:length(vv)));
          dataCount = dataCount + length(vv);
          
        end

        % Now we use full matrices, we can use matmul (opt, R2018a)
        if ip.Results.calcE
          Etheta = Etheta + Nn * ...
            ( YphiExpf*((1i)^(nn+1)*aidx + (-1i)^(nn+1)*pidx) ...
            + YthetaExpf*((1i)^nn*bidx + (-1i)^nn*qidx) );
          Ephi = Ephi + Nn * ...
            (-YthetaExpf*((1i)^(nn+1)*aidx + (-1i)^(nn+1)*pidx) ...
            + YphiExpf*((1i)^nn*bidx + (-1i)^nn*qidx) );
        end
        
        if ip.Results.calcH
          Htheta = Etheta + Nn * ...
            ( YphiExpf*((1i)^(nn+1)*bidx + (-1i)^(nn+1)*qidx) ...
            + YthetaExpf*((1i)^nn*aidx + (-1i)^nn*pidx) );
          Hphi = Ephi + Nn * ...
            (-YthetaExpf*((1i)^(nn+1)*bidx + (-1i)^(nn+1)*qidx) ...
            + YphiExpf*((1i)^nn*aidx + (-1i)^nn*pidx) );
        end
      end

      E=[zeros(size(Etheta)),Etheta,Ephi].';
      H=[zeros(size(Htheta)),Htheta,Hphi].';

      % SI-ify units of H
      H = H * -1i;
    end

    function [E, H, data] = emFieldRtp(beam, rtp, varargin)
      p = inputParser;
  
      p.parse(varargin{:});

      % Scale the locations by the wave number (unitless coordinates)
      rtp(1, :) = rtp(1, :) * abs(beam.k_medium);

      % Get the indices required for the calculation
      cidx = p.Results.cidx;
      if isempty(cidx)
        cidx = find(abs(beam.a)|abs(beam.b));
      end
      [n,m]=ott.utils.combined_index(cidx);
      nm = [ n; m ];

      ci = ott.utils.combined_index(n, m);
      [a, b] = beam.getCoefficients(ci);

      % Calculate the fields
      [E, H, data] = ott.utils.emField(rtp.', beam.basis, nm, [a; b], ...
          'saveData', p.Results.saveData, ...
          'data', p.Results.data, ...
          'calcE', p.Results.calcE, 'calcH', p.Results.calcH);

      % Convert from spherical to Cartesian coordinates
      switch p.Results.coord
        case 'cartesian'
          E = ott.utils.rtpv2xyzv(E,rtp.');
          E(isnan(E)) = 0;
          H = ott.utils.rtpv2xyzv(H,rtp.');
          H(isnan(H)) = 0;

        case 'spherical'
          % Nothing to do

        otherwise
          error('Unknown coordinate system for output');
      end
      
      % Make the matrices 3xN
      E = E.';
      H = H.';
    end

    function [E, H, data] = emFieldXyz(beam, xyz, varargin)
      rtp = ott.utils.xyz2rtp(xyz.');
      [E, H, data] = beam.emFieldRtp(rtp.', varargin{:});
    end
    
    function varargout = visualiseFarfieldSlice(beam, phi, varargin)
      p = inputParser;
      p.addParameter('field', 'irradiance');
      p.addParameter('normalise', false, @islogical);
      p.addParameter('ntheta', 100, @isnumeric);
      p.addParameter('showVisualisation', nargout == 0, @islogical);
      p.parse(varargin{:});
      
      ptheta = linspace(0, 2*pi, p.Results.ntheta);
      
      % TODO: Other field types

      assert(isnumeric(phi), 'phi must be numeric');
      assert(isscalar(phi), ...
        'phi must be scalar in this version, will change in OTT 2.0');

      % Calculate electric field
      [E, ~] = beam.farfield(ptheta, phi);
      
      % Calculate desired field
      [rtp{1:3}] = ott.utils.matchsize(0, ptheta(:), phi);
      I = beam.GetVisualisationData(p.Results.field, [], [rtp{1}, rtp{2}, rtp{3}], [], E.');
%       I = sum(abs(E).^2, 1);
      
      if p.Results.normalise
        I = I ./ max(abs(I(:)));
      end

      % Setup outputs
      if nargout == 2
        varargout{1} = ptheta;
        varargout{2} = I;
      end
      
      % Display visualisation
      if p.Results.showVisualisation
        polarplot(ptheta, I);
      end
      
    end
    
    function visualiseFarfieldSphere(beam, varargin)
      % Generate a spherical surface visualisation of the far-field
      %
      % beam.visualiseFarfieldSphere(phi)
      %
      % Optional named arguments:
      %   npts      num   Number of points to use for sphere surface
      %   normalise bool  If intensity values should be normalised to 1
      %   type      str   Type of visualisation to produce.
      %       sphere    (default) draw a sphere with intensity as color
      %       3dpolar   scale the radius by the intensity
      
      p = inputParser;
      p.addParameter('field', 'irradiance');
      p.addParameter('npts', 100);
      p.addParameter('normalise', false);
      p.addParameter('type', 'sphere');
      p.parse(varargin{:});
      
      % build grid:
      [x,y,z]=sphere(p.Results.npts);

      % generate angular points for farfield:
      [~,theta,phi]=ott.utils.xyz2rtp(x,y,z);

      % find far-field in theta, phi:
      [E,~]=beam.farfield(theta(:),phi(:));
      
      % Calculate the requested field
      dataout = beam.GetVisualisationData(p.Results.field, [], ...
        [theta(:), phi(:), ones(size(phi(:)))], [], E.');

      % Reshape to match the input
      I=reshape(dataout,size(x));
      
      if p.Results.normalise
        I = I ./ max(abs(I(:)));
      end
      
      switch p.Results.type
        case 'sphere'
          surf(x,y,z,I,'facecolor','interp','edgecolor','none');
        case '3dpolar'
          surf(abs(I).*x,abs(I).*y,abs(I).*z,I,...
              'facecolor','interp','edgecolor','none');
        otherwise
          error('Unknown visualisation type');
      end

      zlabel('Z');
      xlabel('X');
      ylabel('Y');
      view(50, 20);
      axis equal;
    end

    function varargout = visualiseFarfield(beam, varargin)
      p = inputParser;
      p.addParameter('size', [80, 80]);
      p.addParameter('direction', 'pos');
      p.addParameter('field', 'irradiance');
      p.addParameter('mapping', 'sin');
      p.addParameter('range', [1, 1]);
      p.addParameter('saveData', nargout == 2);
      p.addParameter('data', []);
      p.addParameter('thetaMax', []);
      p.addParameter('showVisualisation', nargout == 0);
      p.parse(varargin{:});
      
      % If direction is a vector, rotate to that direction
      if ~ischar(p.Results.direction)
        dir = p.Results.direction;
        if numel(dir) == 2
          rbeam = beam.rotateYz(dir(1), dir(2));
        elseif all(size(dir) == [3, 3])
          rbeam = beam.rotate(dir);
        else
          error('OTT:BSC:visualiseFarfield:bad_direction', ...
            'Direction must be char array or 2 element vector or 3x3 matrix');
        end
        
        [varargout{1:nargout}] = rbeam.visualiseFarfield(...
          'size', p.Results.size, 'direction', 'pos', ...
          'field', p.Results.field, 'mapping', p.Results.mapping, ...
          'range', p.Results.range, 'saveData', p.Results.saveData, ...
          'data', p.Results.data, 'thetaMax', p.Results.thetaMax, ...
          'showVisualisation', p.Results.showVisualisation);
        return;  % All done
      end

      % Calculate image locations
      xrange = linspace(-1, 1, p.Results.size(1))*p.Results.range(1);
      yrange = linspace(-1, 1, p.Results.size(2))*p.Results.range(2);
      [xx, yy] = meshgrid(xrange, yrange);

      % Calculate spherical coordinates for pixels
      phi = atan2(yy, xx);
      rr = sqrt(xx.^2 + yy.^2);
      switch p.Results.mapping
        case 'sin'
          theta = asin(rr);
        case 'tan'
          theta = atan(rr);
        case 'theta'
          theta = rr;
        otherwise
          error('Unknown mapping argument value, must be sin, tan or theta');
      end
      
      % Only include points within NA range
      thetaMax = Inf;
      if ~isempty(p.Results.thetaMax)
        thetaMax = p.Results.thetaMax;
      end

      % Determine if the points need calculating
      pinside = imag(theta) == 0 & theta < thetaMax;
      iphi = phi(pinside);
      itheta = theta(pinside);
      
      if strcmpi(p.Results.direction, 'neg')
        itheta = pi - itheta;
      elseif ~strcmpi(p.Results.direction, 'pos')
        error('Direction must be \'pos\' or \'neg\'');
      end

      % Calculate the electric field in the farfield
      [ioutputE, ~, data] = beam.farfield(itheta(:), iphi(:), ...
        'saveData', p.Results.saveData, 'data', p.Results.data, ...
        'calcE', true, 'calcH', false);
      
      % Generate the requested field
      dataout = beam.GetVisualisationData(p.Results.field, [], ...
        [itheta, iphi, ones(size(iphi))], [], ioutputE.');

      % Pack the result into the images
      imout = zeros(p.Results.size);
      imout(pinside) = dataout;

      % Display the visualisation
      if p.Results.showVisualisation
        
        % Check the field is real
        if ~isreal(imout)
          error(['Unsupported field type for visualisation: ' p.Results.field]);
        end
        
        imagesc(xrange, yrange, imout);
        caxis([min(dataout), max(dataout)]);
        xlabel('X');
        ylabel('Y');
        axis image;
      end
      
      % Handle outputs
      if nargout == 1
        varargout{1} = imout;
      elseif nargout == 2
        varargout{1} = imout;
        varargout{2} = data;
      end
    end
'''

'''
    function varargout = visualise(beam, varargin)
      p = inputParser;
      p.addParameter('field', 'irradiance');
      p.addParameter('size', [ 80, 80 ]);
      p.addParameter('axis', 'z');
      p.addParameter('offset', 0.0);
      p.addParameter('range', ...
          [1,1]*ott.utils.nmax2ka(beam.Nmax)/abs(beam.k_medium));
      p.addParameter('saveData', nargout == 2);
      p.addParameter('data', []);
      p.addParameter('mask', []);
      p.addParameter('showVisualisation', nargout == 0);
      p.addParameter('combine', []);
      p.addParameter('cidx', []);   % Used with savedata
      p.parse(varargin{:});

      if iscell(p.Results.range)
        xrange = p.Results.range{1};
        yrange = p.Results.range{2};
        sz = [length(yrange), length(xrange)];
      elseif length(p.Results.range) == 2
        xrange = linspace(-1, 1, p.Results.size(1))*p.Results.range(1);
        yrange = linspace(-1, 1, p.Results.size(2))*p.Results.range(2);
        sz = p.Results.size;
      elseif length(p.Results.range) == 4
        xrange = linspace(p.Results.range(1), p.Results.range(2), p.Results.size(1));
        yrange = linspace(p.Results.range(3), p.Results.range(4), p.Results.size(2));
        sz = p.Results.size;
      else
        error('ott:Bsc:visualise:range_error', 'Incorrect number of range arguments');
      end
      [xx, yy, zz] = meshgrid(xrange, yrange, p.Results.offset);

      % Generate the xyz grid for the used requested plane
      if ischar(p.Results.axis)
        switch p.Results.axis
          case 'x'
            xyz = [zz(:), yy(:), xx(:)];
            alabels = {'Z', 'Y'};
          case 'y'
            xyz = [yy(:), zz(:), xx(:)];
            alabels = {'Z', 'X'};
          case 'z'
            xyz = [xx(:), yy(:), zz(:)];
            alabels = {'X', 'Y'};
          otherwise
            error('Unknown axis name specified');
        end
      elseif iscell(p.Results.axis)
        dir1 = p.Results.axis{1}(:);
        dir2 = p.Results.axis{2}(:);
        if numel(p.Results.axis) == 3
          dir3 = p.Results.axis{3}(:);
        else
          dir3 = cross(dir1(:), dir2(:));
        end
        
        alabels = {'Direction 1', 'Direction 2'};
        
        xyz = dir1.*xx(:).' + dir2.*yy(:).' + dir3.*zz(:).';
        xyz = xyz.';
      else
        error('axis must be character or cell array');
      end
      
      if strcmpi(p.Results.combine, 'coherent')
        % Reduce the amount of work done in emFieldXyz
        beam = sum(beam);
      end
      
      imout = zeros(sz(2), sz(1), beam.Nbeams);
      
      % If save data is requested, this gets reused for multiple beams
      data = p.Results.data;
      
      for ii = 1:beam.Nbeams

        % Calculate the electric field
        [E, ~, data] = beam.beam(ii).emFieldXyz(xyz.', ...
            'saveData', p.Results.saveData, 'data', data, ...
            'calcE', true', 'calcH', false, 'cidx', p.Results.cidx);

        % Generate the requested field
        dataout = beam.GetVisualisationData(p.Results.field, ...
          xyz, [], E.', []);

        % Reshape the output
        imout(:, :, ii) = reshape(dataout, sz(2), sz(1));
        
      end
      
      if strcmpi(p.Results.combine, 'incoherent')
        imout = sum(imout, 3);
      end

      % Display the visualisation
      if p.Results.showVisualisation && ismatrix(imout)
        
        % Apply the mask
        if ~isempty(p.Results.mask)
          imout(p.Results.mask(xyz.')) = NaN;
        end
        
        imagesc(xrange, yrange, imout, 'AlphaData', ~isnan(imout));
        xlabel(alabels{1});
        ylabel(alabels{2});
        axis image;
      elseif p.Results.showVisualisation
        warning('OTT:Bsc:visualise:too_many_beams', ...
          ['Visualisation not shown for multiple beams\n', ...
           '> Consider combining beams coherently or incoherently']);
      end
      
      % Handle outputs
      if nargout == 1
        varargout{1} = imout;
      elseif nargout == 2
        varargout{1} = imout;
        varargout{2} = data;
      end
    end

    function p = get.power(beam)
      % get.power calculate the power of the beam
      p = full(sum(abs(beam.a).^2 + abs(beam.b).^2));
    end

    function beam = set.power(beam, p)
      % set.power set the beam power
      beam = sqrt(p / beam.power) * beam;

      warning('ott:Bsc:set_power_not_recomended', ...
          ['Changing the beam power in this or previous versions', newline, ...
          'of OTT not recomended, see documentation for more info']);
    end

    function lambda = get.wavelength(beam)
      % Get the beam wavelength
      lambda = 2*pi/beam.k_medium;
    end

    function beam = set.wavelength(beam, lambda)
      % Set the beam wavelength
      beam.k_medium = 2*pi/lambda;
    end

    function speed = get.speed(beam)
      % Get the speed of the beam in medium
      speed = beam.omega / beam.k_medium;
    end

    function beam = set.basis(beam, basis)
      % Set the beam type, checking it is a valid type first
      if ~any(strcmpi(basis, {'incoming', 'outgoing', 'regular'}))
        error('OTT:Bsc:set_basis:invalid_value', 'Invalid beam basis');
      end
      beam.basis = basis;
    end

    function beam = set.type(beam, type)
      % Set the beam type, checking it is a valid type first
      if ~any(strcmpi(type, {'incident', 'scattered', 'total', 'internal'}))
        error('OTT:Bsc:set_type:invalid_value', 'Invalid beam type');
      end
      beam.type = type;
    end

    function nbeams = get.Nbeams(beam)
      % get.beams get the number of beams in this object
      nbeams = size(beam.a, 2);
    end

    function bsc = beam(bsc, idx)
      % BEAM get beams from a beam array object
      %
      % BEAM(idx) idx can be a linear index or a logical array.
      bsc.a = bsc.a(:, idx);
      bsc.b = bsc.b(:, idx);
    end

    function nmax = get.Nmax(beam)
      %get.Nmax calculates Nmax from the current size of the beam coefficients
      nmax = ott.utils.combined_index(size(beam.a, 1));
    end

    function beam = set.Nmax(beam, nmax)
      %set.Nmax resizes the beam vectors
      beam = beam.set_Nmax(nmax);
    end

    function nbeam = shrink_Nmax(beam, varargin)
      % SHRINK_NMAX reduces the size of the beam while preserving power

      p = inputParser;
      p.addParameter('tolerance', 1.0e-6);
      p.parse(varargin{:});

      amagA = full(sum(sum(abs(beam.a).^2)));
      bmagA = full(sum(sum(abs(beam.b).^2)));

      for ii = 1:beam.Nmax

        total_orders = ott.utils.combined_index(ii, ii);
        nbeam = beam;
        nbeam.a = nbeam.a(1:total_orders);
        nbeam.b = nbeam.b(1:total_orders);

        amagB = full(sum(sum(abs(nbeam.a).^2)));
        bmagB = full(sum(sum(abs(nbeam.b).^2)));

        aapparent_error = abs( amagA - amagB )/amagA;
        bapparent_error = abs( bmagA - bmagB )/bmagA;

        if aapparent_error < p.Results.tolerance && ...
            bapparent_error < p.Results.tolerance
          break;
        end
      end
    end

    function beam = set_Nmax(beam, nmax, varargin)
      p = inputParser;
      p.addParameter('tolerance', 1.0e-6);
      p.addParameter('powerloss', 'warn');
      p.parse(varargin{:});

      total_orders = ott.utils.combined_index(nmax, nmax);
      if size(beam.a, 1) > total_orders

        amagA = full(sum(sum(abs(beam.a).^2)));
        bmagA = full(sum(sum(abs(beam.b).^2)));

        beam.a = beam.a(1:total_orders, :);
        beam.b = beam.b(1:total_orders, :);

        amagB = full(sum(sum(abs(beam.a).^2)));
        bmagB = full(sum(sum(abs(beam.b).^2)));

        if ~strcmpi(p.Results.powerloss, 'ignore')

          aapparent_error = abs( amagA - amagB )/amagA;
          bapparent_error = abs( bmagA - bmagB )/bmagA;

          if aapparent_error > p.Results.tolerance || ...
              bapparent_error > p.Results.tolerance
            if strcmpi(p.Results.powerloss, 'warn')
              warning('ott:Bsc:setNmax:truncation', ...
                  ['Apparent errors of ' num2str(aapparent_error) ...
                      ', ' num2str(bapparent_error) ]);
            elseif strcmpi(p.Results.powerloss, 'error')
              error('ott:Bsc:setNmax:truncation', ...
                  ['Apparent errors of ' num2str(aapparent_error) ...
                      ', ' num2str(bapparent_error) ]);
            else
              error('ott:Bsc:setNmax:truncation', ...
                'powerloss should be one of ignore, warn or error');
            end
          end
        end
      elseif size(beam.a, 1) < total_orders
        [arow_index,acol_index,aa] = find(beam.a);
        [brow_index,bcol_index,ba] = find(beam.b);
        beam.a = sparse(arow_index,acol_index,aa,total_orders,beam.Nbeams);
        beam.b = sparse(brow_index,bcol_index,ba,total_orders,beam.Nbeams);
      end
    end

    function beam = translate(beam, A, B)
      % TRANSLATE apply a translation using given translation matrices.
      %
      % TRANSLATE(A, B) applies the translation given by A, B.
      beam = [ A B ; B A ] * beam;
    end

    function [beam, A, B] = translateZ(beam, varargin)
      p = inputParser;
      p.addOptional('z', []);
      p.addParameter('Nmax', beam.Nmax);
      p.parse(varargin{:});

      if nargout ~= 1 && numel(p.Results.z) > 1
        error('Multiple output with multiple translations not supported');
      end

      if ~isempty(p.Results.z)
        z = p.Results.z;

        % Add a warning when the beam is translated outside nmax2ka(Nmax) 
        % The first time may be OK, the second time does not have enough
        % information.
        if beam.dz > ott.utils.nmax2ka(beam.Nmax)/beam.k_medium
          warning('ott:Bsc:translateZ:outside_nmax', ...
              'Repeated translation of beam outside Nmax region');
        end
        beam.dz = beam.dz + abs(z);

        % Convert to beam units
        z = z * beam.k_medium / 2 / pi;

        ibeam = beam;
        beam = ott.Bsc();

        for ii = 1:numel(z)
          [A, B] = ibeam.translateZ_type_helper(z(ii), [p.Results.Nmax, ibeam.Nmax]);
          beam = beam.append(ibeam.translate(A, B));
          beam.basis = 'regular';
        end
      else
        error('Wrong number of arguments');
      end

      % Pack the rotated matricies into a single ABBA object
      if nargout == 2
        A = [ A B ; B A ];
      end
    end

    function varargout = translateXyz(beam, varargin)
      p = inputParser;
      p.addOptional('opt1', []);    % xyz or Az
      p.addOptional('opt2', []);    % [] or Bz
      p.addOptional('opt3', []);    % [] or D
      p.addParameter('Nmax', beam.Nmax);
      p.parse(varargin{:});

      if ~isempty(p.Results.opt1) && isempty(p.Results.opt2) ...
          && isempty(p.Results.opt3)
        xyz = p.Results.opt1;
        rtp = ott.utils.xyz2rtp(xyz.').';
        [varargout{1:nargout}] = beam.translateRtp(rtp, ...
            'Nmax', p.Results.Nmax);
      else
        [varargout{1:nargout}] = beam.translateRtp(varargin{:});
      end
    end

    function [beam, A, B, D] = translateRtp(beam, varargin)
      p = inputParser;
      p.addOptional('opt1', []);    % rtp or Az
      p.addOptional('opt2', []);    % [] or Bz
      p.addOptional('opt3', []);    % [] or D
      p.addParameter('Nmax', beam.Nmax);
      p.parse(varargin{:});

      % Convert Nmax to a single number
      if numel(p.Results.Nmax) == 1
        oNmax = p.Results.Nmax;
      elseif numel(p.Results.Nmax) == 2
        oNmax = p.Results.Nmax(2);
      else
        error('Nmax must be 2 element vector or scalar');
      end

      % Handle input arguments
      if ~isempty(p.Results.opt1) && isempty(p.Results.opt2) ...
          && isempty(p.Results.opt3)

        % Assume first argument is rtp coordinates
        r = p.Results.opt1(1, :);
        theta = p.Results.opt1(2, :);
        phi = p.Results.opt1(3, :);

      elseif ~isempty(p.Results.opt1) && ~isempty(p.Results.opt2) ...
          && ~isempty(p.Results.opt3)

        % Rotation/translation is already computed, apply it
        A = p.Results.opt1;
        B = p.Results.opt2;
        D = p.Results.opt3;
        beam = beam.rotate('wigner', D);
        beam = beam.translate(A, B);
        beam = beam.rotate('wigner', D');
        return;
      else
        error('Not enough input arguments');
      end

      if numel(r) ~= 1 && nargout ~= 1
        error('Multiple output with multiple translations not supported');
      end

      % Only do the rotation if we need it
      if any((theta ~= 0 & abs(theta) ~= pi) | phi ~= 0)

        ibeam = beam;
        beam = ott.Bsc();

        for ii = 1:numel(r)
          [tbeam, D] = ibeam.rotateYz(theta(ii), phi(ii), ...
              'Nmax', max(oNmax, ibeam.Nmax));
          [tbeam, A, B] = tbeam.translateZ(r(ii), 'Nmax', oNmax);
          beam = beam.append(tbeam.rotate('wigner', D'));
        end
      else
        dnmax = max(oNmax, beam.Nmax);
        D = speye(ott.utils.combined_index(dnmax, dnmax));

        % Replace rotations by 180 with negative translations
        idx = abs(theta) == pi;
        r(idx) = -r(idx);

        if numel(r) == 1
          [beam, A, B] = beam.translateZ(r, 'Nmax', oNmax);
        else
          beam = beam.translateZ(r, 'Nmax', oNmax);
        end
      end

      % Rotate the translation matricies
      if nargout == 3 || nargout == 2

        % The beam might change size, so readjust D to match
        sz = size(A, 1);
        D2 = D(1:sz, 1:sz);

        A = D2' * A * D;
        B = D2' * B * D;

        % Pack the rotated matricies into a single ABBA object
        if nargout == 2
          A = [ A B ; B A ];
        end
      elseif nargout ~= 4 && nargout ~= 1
        error('Insufficient number of output arguments');
      end
    end

    function [beam, D] = rotate(beam, varargin)
      p = inputParser;
      p.addOptional('R', []);
      p.addParameter('Nmax', beam.Nmax);
      p.addParameter('wigner', []);
      p.parse(varargin{:});

      if ~isempty(p.Results.R) && isempty(p.Results.wigner)

        R = p.Results.R;

        % If no rotation, don't calculate wigner rotation matrix
        if sum(sum((eye(3) - R).^2)) < 1e-6
          D = eye(ott.utils.combined_index(p.Results.Nmax, p.Results.Nmax));
          return;
        end

        D = ott.utils.wigner_rotation_matrix(...
            max(beam.Nmax, p.Results.Nmax), R);
        beam = beam.rotate('wigner', D);

      elseif ~isempty(p.Results.wigner) && isempty(p.Results.R)
        
        if iscell(p.Results.wigner)
          
          ibeam = beam;
          beam = ott.Bsc();

          for ii = 1:numel(p.Results.wigner)
            sz = size(ibeam.a, 1);
            D2 = p.Results.wigner{ii}(1:sz, 1:sz);
            beam = beam.append(D2 * ibeam);
          end
        
        else

          sz = size(beam.a, 1);
          D2 = p.Results.wigner(1:sz, 1:sz);
          beam = D2 * beam;
          
        end

      else
        error('One of wigner or R must be specified');
      end
    end

    function [beam, D] = rotateX(beam, angle, varargin)
      %ROTATEX rotates the beam about the x-axis an angle in radians
      import ott.utils.*;
      [beam, D] = beam.rotate(rotx(angle*180/pi), varargin{:});
    end

    function [beam, D] = rotateY(beam, angle, varargin)
      %ROTATEX rotates the beam about the y-axis an angle in radians
      import ott.utils.*;
      [beam, D] = beam.rotate(roty(angle*180/pi), varargin{:});
    end

    function [beam, D] = rotateZ(beam, angle, varargin)
      %ROTATEX rotates the beam about the z-axis an angle in radians
      import ott.utils.*;
      [beam, D] = beam.rotate(rotz(angle*180/pi), varargin{:});
    end

    function [beam, D] = rotateXy(beam, anglex, angley, varargin)
      %ROTATEX rotates the beam about the x then y axes
      import ott.utils.*;
      [beam, D] = beam.rotate(roty(angley*180/pi)*rotx(anglex*180/pi), ...
          varargin{:});
    end

    function [beam, D] = rotateXz(beam, anglex, anglez, varargin)
      %ROTATEX rotates the beam about the x then z axes
      import ott.utils.*;
      [beam, D] = beam.rotate(rotz(anglez*180/pi)*rotx(anglex*180/pi), ...
          varargin{:});
    end

    function [beam, D] = rotateYz(beam, angley, anglez, varargin)
      %ROTATEX rotates the beam about the y then z axes
      import ott.utils.*;
      [beam, D] = beam.rotate(rotz(anglez*180/pi)*roty(angley*180/pi), ...
          varargin{:});
    end

    function [beam, D] = rotateXyz(beam, anglex, angley, anglez, varargin)
      import ott.utils.*;
      [beam, D] = beam.rotate(rotz(anglez*180/pi)* ...
          roty(angley*180/pi)*rotx(anglex*180/pi), varargin{:});
    end

    function beam = totalField(beam, ibeam)
      switch beam.type
        case 'total'
          % Nothing to do
          
        case 'scattered'
          beam = 2*beam + ibeam;
          beam.type = 'total';
          
        case 'internal'
          error('Cannot convert from internal to total field');
          
        case 'incident'
          error('Cannot convert from incident to total field');
          
        otherwise
          error('Unknown beam type');
      end
    end

    def scattered_field(self, ibeam)
      % Calculate the scattered field representation of the beam
      %
      % scattered_beam = beam.totalField(incident_beam)
      if self.type == 'total':
        
        case 'total'
                  beam = 0.5*(beam - ibeam);
          beam.type = 'scattered';
      if self.type == 'total':
          
        case 'scattered'
          % Nothing to do
      if self.type == 'total':
          
        case 'internal'
          error('Cannot convert from internal to scattered field');
      if self.type == 'total':
          
        case 'incident'
          error('Cannot convert from incident to total field');
      if self.type == 'total':
          
        otherwise
          error('Unknown beam type');

    function [a, b] = getCoefficients(beam, ci)
      if nargin == 1
        ci = 1:size(beam.a, 1);
      end

      a = beam.a(ci, :);
      b = beam.b(ci, :);

      if nargout == 1
        a = [a; b];
      end
    end

    function [n, m] = getModeIndices(beam)
      %GETMODEINDICES gets the mode indices
      [n, m] = ott.utils.combined_index([1:size(beam.a, 1)].');
      if nargout == 1
        n = [n; m];
      end
    end
'''




'''
    def divide(self, scalar):
        self.a = self.a/scalar
        self.b = self.b/scalar
    
    function [moment, int, data] = intensityMoment(beam, varargin)
      p = inputParser;
      p.addParameter('thetaRange', [0, pi]);
      p.addParameter('saveData', false);
      p.addParameter('ntheta', 100);
      p.addParameter('nphi', 100);
      p.addParameter('data', []);
      p.parse(varargin{:});

      % Regular beams have trivial solution
      if strcmpi(beam.basis, 'regular')
        moment = [0;0;0];
        int = 0;
        data = p.Results.data;
        warning('Regular wavefunctions go to zero in far-field');
        return;
      end

      % Setup the angular grid
      [theta, phi] = ott.utils.angulargrid(p.Results.ntheta, p.Results.nphi);
      dtheta = theta(2) - theta(1);
      dphi = phi(p.Results.ntheta+1) - phi(1);

      % Truncate the theta range
      keep = theta > p.Results.thetaRange(1) & theta < p.Results.thetaRange(2);
      theta = theta(keep);
      phi = phi(keep);

      uxyz = ott.utils.rtp2xyz([ones(size(theta)), theta, phi]).';
      
      % So integrals match sign convention used in ott.forcetorque
      uxyz(3, :) = -uxyz(3, :);

      % Calculate E-field in far-field
      [E, ~, data] = beam.farfield(theta, phi, ...
          'saveData', p.Results.saveData, ...
          'data', p.Results.data);
'''


'''
      % Calculate the irradiance
      Eirr = sum(abs(E).^2, 1);
      int = sum(Eirr .* sin(theta.') .* dtheta .* dphi, 2);

      % Calculate moment in Cartesian coordinates
      Eirr_xyz = uxyz .* Eirr;
      moment = sum(Eirr_xyz .* sin(theta.') .* dtheta .* dphi, 2);
    end

    function [sbeam, beam] = scatter(beam, tmatrix, varargin)
      p = inputParser;
      p.addParameter('position', []);
      p.addParameter('rotation', []);
      p.parse(varargin{:});

      % Determine the maximum tmatrix.Nmax(2) and check type
      maxNmax1 = 0;
      maxNmax2 = 0;
      tType = tmatrix(1).type;
      for ii = 1:numel(tmatrix)
        maxNmax1 = max(maxNmax1, tmatrix(ii).Nmax(1));
        maxNmax2 = max(maxNmax2, tmatrix(ii).Nmax(2));
        if ~strcmpi(tmatrix(ii).type, tType)
          error('T-matrices must be same type');
        end
      end

      % If the T is scattered, we can save time by throwing away columns
      % Only works when we don't grow the beam in translation
      if strcmpi(tmatrix(1).type, 'scattered') ...
          && isempty(p.Results.position)
        maxNmax2 = min(maxNmax2, beam.Nmax);
      end

      % Ensure all T-matrices are the same size
      for ii = 1:numel(tmatrix)
        tmatrix(ii).Nmax = [maxNmax1, maxNmax2];
      end

      % Apply translation to the beam
      if ~isempty(p.Results.position)

        % Requires scattered beam, convert if needed
        if ~strcmpi(tmatrix(1).type, 'scattered')
          maxNmax2 = min(maxNmax2, beam.Nmax);
          for ii = 1:numel(tmatrix)
            tmatrix(ii).type = 'scattered';
            tmatrix(ii).Nmax = [maxNmax1, maxNmax2];
          end
        end

'''

'''
        % Apply translation
        % We need Nmax+1 terms for the force calculation
        beam = beam.translateXyz(p.Results.position, 'Nmax', maxNmax2+1);
      end

      % Apply rotation to the beam
      rbeam = beam;
      if ~isempty(p.Results.rotation)
        [rbeam, D] = rbeam.rotate(p.Results.rotation, ...
            'Nmax', maxNmax1);
      end

      % Ensure the Nmax for the inner dimension matches
      if strcmpi(tmatrix(1).type, 'scattered')
        % T-matrix is already done
        rbeam = rbeam.set_Nmax(maxNmax2, 'powerloss', 'ignore');
      else
        for ii = 1:numel(tmatrix)
          tmatrix(ii) = tmatrix(ii).set_Nmax([maxNmax1, rbeam.Nmax], ...
              'powerloss', 'ignore');
        end
        if ~strcmpi(tmatrix(1).type, 'internal')
          ott.warning('ott:Bsc:scatter', ...
              'It may be more optimal to use a scattered T-matrix');
        end
      end

      % Calculate the resulting beams
      sbeam = ott.Bsc();
      for ii = 1:numel(tmatrix)
        sbeam = sbeam.append(tmatrix(ii).data * rbeam);
      end

      % Apply the inverse rotation
      if ~isempty(p.Results.rotation)

        % This seems to take a long time
        %sbeam = sbeam.rotate('wigner', D');

        sbeam = sbeam.rotate(inv(p.Results.rotation));
      end

      % Assign a type to the resulting beam
      switch tmatrix(1).type
        case 'total'
          sbeam.type = 'total';
          sbeam.basis = 'outgoing';
          
        case 'scattered'
          sbeam.type = 'scattered';
          sbeam.basis = 'outgoing';
          
        case 'internal'
          sbeam.type = 'internal';
          sbeam.basis = 'regular';
        
          % Wavelength has changed, update it
          sbeam.k_medium = tmatrix(1).k_particle;
          
        otherwise
          error('Unrecognized T-matrix type');
      end
    end

    function beam = mtimes(a,b)
      if isa(a, 'ott.Bsc')
        beam = a;
        beam.a = beam.a * b;
        beam.b = beam.b * b;
      else
        beam = b;
        if size(a, 2) == 2*size(beam.a, 1)
          ab = a * [beam.a; beam.b];
          beam.a = ab(1:size(ab, 1)/2, :);
          beam.b = ab(1+size(ab, 1)/2:end, :);
        else
          beam.a = a * beam.a;
          beam.b = a * beam.b;
        end
      end
    end

    function beam = plus(beam1, beam2)
      %PLUS add two beams together

      if beam1.Nmax > beam2.Nmax
        beam2.Nmax = beam1.Nmax;
      elseif beam2.Nmax > beam1.Nmax
        beam1.Nmax = beam2.Nmax;
      end

      beam = beam1;
      beam.a = beam.a + beam2.a;
      beam.b = beam.b + beam2.b;
    end

    function beam = minus(beam1, beam2)
      %MINUS subtract two beams

      if beam1.Nmax > beam2.Nmax
        beam2.Nmax = beam1.Nmax;
      elseif beam2.Nmax > beam1.Nmax
        beam1.Nmax = beam2.Nmax;
      end

      beam = beam1;
      beam.a = beam.a - beam2.a;
      beam.b = beam.b - beam2.b;
    end
    
    function beam = sum(beamin, dim)
      if numel(beamin) > 1
        % beam is an array
        
        % Handle default value for dimension
        if nargin < 2
          if isvector(beamin)
            dim = find(size(beamin) > 1, 1);
          elseif ismatrix(beamin) == 2
            dim = 2;
          else
            dim = find(size(beamin) > 1, 1);
          end
        end
        
        % Select the first row in our dimension
        subs = [repmat({':'},1,dim-1), 1, ...
          repmat({':'},1,ndims(beamin)-dim)];
        S = struct('type', '()', 'subs', {subs});
        beam = subsref(beamin, S);
        
        % Add each beam
        for ii = 2:size(beamin, dim)
        subs = [repmat({':'},1,dim-1), ii, ...
          repmat({':'},1,ndims(beamin)-dim)];
          S = struct('type', '()', 'subs', {subs});
          beam = beam + subsref(beamin, S);
        end
        
      else
        % Beam union
        beam = beamin;
        beam.a = sum(beam.a, 2);
        beam.b = sum(beam.b, 2);
      end
    end
    
    def clear_dz(self):
        self.dz = 0
'''