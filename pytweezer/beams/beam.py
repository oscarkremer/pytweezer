import warnings
import numpy as np
from pytweezer.utils import combined_index, translate_z, xyz2rtp
from .translation import translate_xyz
from scipy.sparse import csr_matrix
from copy import copy
import time

class Beam:
    def __init__(self, a, b, basis, beam_type, k_m=2*np.pi, omega=2*np.pi, dz=0):
        pass
#        self.dz = dz
#        self.k_medium = k_m
#        self.omega = omega
#
#        assert a.shape==b.shape, 'Mismatch between the shapes of the coefficients'      
#        if isinstance(a, np.ndarray):
#            if len(a.shape) == 1:
#                a = a.reshape((a.shape[0], 1))
#        if isinstance(b, np.ndarray):
#            if len(b.shape) == 1:
#                b = b.reshape((b.shape[0], 1))
#        assert a.shape[0] >= 3, 'Number of multiploe must be at least 3'
#        assert np.sqrt(a.shape[0]+1) == np.floor(np.sqrt(a.shape[0]+1)), 'Number of multipoles must be 3, 8, 15, 24, ...'
#        self.a = a
#        self.b = b
#        self.basis = basis
#        self.beam_type = beam_type

    def compute_k_medium(self, index_m, lambda_0):
        if index_m < 0 or lambda_0 < 0:
            raise ValueError('Wavelength and refractive index cannot be negative.')
        else:
            self.k_m = index_m*2*np.pi/lambda_0
            self.index_m = index_m

    def make_beam_vector(self, a, b, n, m, n_max=None):
        if not n_max:
            if not n.size:
                self.n_max = 0
            else:
                self.n_max = n.max()
            n_max = self.n_max
        total_orders = combined_index(n_max, n_max)
        ci = combined_index(n, m)
        n_beams = a.shape[1] if len(a.shape) > 1 else 1
        ci, c_in_beams = np.meshgrid(ci, np.arange(1, n_beams+1, 1), indexing='ij')
        ci = ci.reshape((-1))
        c_in_beams = c_in_beams.reshape((-1))
        a = csr_matrix((a, (ci-1, c_in_beams-1)), shape=(total_orders, n_beams)).toarray()
        b = csr_matrix((b, (ci-1, c_in_beams-1)), shape=(total_orders, n_beams)).toarray() 
        
        n, m = combined_index(np.arange(1, n_max**2+2*n_max+1, 1))
        n = n.T
        m = m.T
        return a, b, n, m

    def translate_z_type_helper(self, z, n_max):
        if self.basis == 'incoming':
            translation_type = 'sbesselh2';
        elif self.basis == 'outgoing':
            translation_type = 'sbesselh1';
        elif self.basis == 'regular':
            translation_type = 'sbesselj'
        A, B = translate_z(n_max, z, function_type=translation_type)
        return A, B

    def shrink_nmax(self, tolerance=1e-6):
        amagA = np.power(np.abs(self.a), 2).sum()
        bmagA = np.power(np.abs(self.b), 2).sum()
        for i in range(1, self.n_max+1):
            total_orders = combined_index(i, i)
            new_beam = copy(self)
            new_beam.a = new_beam.a[:total_orders]
            new_beam.b = new_beam.b[:total_orders]
            amagB = np.power(np.abs(new_beam.a), 2).sum()
            bmagB = np.power(np.abs(new_beam.b), 2).sum()
            a_error = abs( amagA - amagB )/amagA
            b_error = abs( bmagA - bmagB )/bmagA
            if a_error < tolerance and b_error < tolerance:
                self.a = new_beam.a
                self.b = new_beam.b 
                self.n_max = i
                break

    def scatter(self, t_matrix, position=np.array([[],[],[]]), rotation=np.array([[],[],[]])):
        max_n_max1 = 0
        max_n_max2 = 0
        t_type = t_matrix.type
        if isinstance(t_matrix, np.ndarray):
            for t in t_matrix:
                maxNmax1 = max(maxNmax1, t.Nmax[0])
                maxNmax2 = max(maxNmax2, t.Nmax[1])
            # TODO: implement method for when t_matrix is an array of matrices
        else:

            max_n_max1 = max(max_n_max1, t_matrix.get_n_max()[0])
            max_n_max2  = max(max_n_max2, t_matrix.get_n_max()[1])
            if t_matrix.type == 'scattered' and not position.size:
                max_n_max2 = min(max_n_max2, self.n_max)
            #make set nmax method for t_matrix
            t_matrix.set_n_max(np.array([max_n_max1, max_n_max2]))
            if position.size:
                if t_matrix.type != 'scattered':
                    max_n_max2 = min(max_n_max2, self.n_max)
                    t_matrix.set_type('scattered')
                    t_matrix.set_n_max(np.array([max_n_max1, max_n_max2]))
                #checked until here

                beam = translate_xyz(copy(self), position, n_max=max_n_max2+1)
            start = time.time()

            r_beam = copy(beam)
            if rotation.size:                            
                r_beam, D = r_beam.rotate(rotation, n_max=max_n_max1)
            if t_matrix.type == 'scattered':
                r_beam.set_n_max(max_n_max2, power_loss='ignore')
            else:
                t_matrix.set_n_max(np.array([max_n_max1, r_beam.n_max]), power_loss='ignore')
                if t_matrix.type == 'internal':
                    warnings.warn('ott:Bsc:scatter: It may be more optimal to use a scattered T-matrix');
            s_beam = t_matrix.T*r_beam
            if rotation.size:
                s_beam = s_beam.rotate(np.linalg.inv(p.Results.rotation))
            elif t_matrix.type == 'total':
                s_beam.type = 'total'
                s_beam.basis = 'outgoing'
            elif t_matrix.type == 'scattered':
                s_beam.set_type('scattered')
                s_beam.set_basis('outgoing') 
            elif t_matrix.type == 'internal':
                s_beam.set_type('internal')# = 'internal'
                s_beam.set_basis('regular')
                s_beam.k_m = tmatrix.k_p
                return NotImplemented
            else:
                raise ValueError('Unrecognized T-matrix type')
            end = time.time()
            return s_beam, beam

    def set_basis(self, basis):
        if isinstance(basis, str):
            if basis in ['incoming', 'outgoing', 'regular']:
                self._basis_ = basis
            else:
                raise ValueError('Basis must be incoming, outgoing or regular')
        else:
            raise TypeError('Basis must be string type')

    def set_type(self, beam_type):
        if isinstance(beam_type, str):
            if beam_type in ['incident', 'scattered', 'total', 'internal']:
                self._type_ = beam_type
            else:
                raise ValueError('Beam type must be incident, scattered, total or internal')
        else:
            raise TypeError('Basis must be string type')

    def total_field(self, i_beam):
        if self._type_ == 'total':
            pass
        elif self._type_ == 'scattered':    
            beam = copy(self)   
            beam = 2*beam + i_beam
            beam.set_type('total')
            return beam
        elif self._type_ == 'internal':
            raise ValueError('Cannot convert from internal to total field')
        elif self._type_ == 'incident':
            raise ValueError('Cannot convert from incident to total field')
        else:
            raise ValueError('Unknown beam type')

    def get_basis(self):
        return self._basis_

    def get_n_max(self):
        return int(combined_index(self.a.shape[0])[0])
        
    def get_n_beams(self):
        return self.a.shape[1]
    
    def set_n_max(self, n_max, tolerance=1e-6, power_loss='warn'):
        total_orders = int(combined_index(n_max, n_max))
        if self.a.shape[0] > total_orders:
            amagA = np.power(np.abs(self.a),2).sum()
            bmagA = np.power(np.abs(self.b),2).sum()
            amagB = np.power(np.abs(self.a[:total_orders, :]),2).sum()
            bmagB = np.power(np.abs(self.b[:total_orders, :]),2).sum()
            if not power_loss == 'ignore':
                aapparent_error = abs( amagA - amagB )/amagA
                bapparent_error = abs( bmagA - bmagB )/bmagA
                if aapparent_error > tolerance or bapparent_error > tolerance:
                    if power_loss == 'warn':
                        warnings.warn('ott:Bsc:setNmax:truncation')
                    elif power_loss=='error':
                        raise ValueError('Truncation error during order reduction')
                    else:
                        raise ValueError('Power loss variable should be: ignore, warn or error')
                else:
                    self.a = self.a[:total_orders, :]
                    self.b = self.b[:total_orders, :]
            else:
                self.a = self.a[:total_orders, :]
                self.b = self.b[:total_orders, :]
        elif self.a.shape[0] < total_orders:
            rows, columns = np.where(self.a)
            elements = self.a[(rows, columns)]
            self.a = csr_matrix((elements, (rows, columns)), shape=(total_orders, self._n_beams_)).toarray()
            rows, columns = np.where(self.b)
            elements = self.b[(rows, columns)]
            self.b = csr_matrix((elements, (rows, columns)), shape=(total_orders, self._n_beams_)).toarray() 
        

    def get_coefficients(self):
        return self.a, self.b

    def get_mode_indices(self):
        n, m = combined_index(np.arange(1, self.a.shape[0]+1).T)
        return n.astype(int), m.astype(int)

    def append(self, other):
        if self._n_beams_ == 1:
            # TODO: Have to implement this section
            self.a = other.a
            self.b = other.b
        else:
            self.n_max = max(self.n_max, other.n_max)
            other.n_max = self.n_max
            self.a = np.concatenate([self.a, other.a])
            self.b = np.concatenate([self.b, other.b])

    def __add__(self, beam2):
        if self.n_max > beam2.n_max:
            beam2.set_n_max(self.n_max)
        elif beam2.n_max > self.n_max:
            self.set_n_max(beam2.n_max)
        else:
            pass
        self.a = self.a + beam2.a
        self.b = self.b + beam2.b
        return self

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            self.a = other*self.a
            self.b = other*self.b
            return self
        elif isinstance(other, np.ndarray):
            if other.shape[1] == 2*self.a.shape[0]:
                ab = np.matmul(other, np.concatenate([self.a, self.b]))
                self.a = ab[:int(ab.shape[0]/2),:]
                self.b = ab[int(ab.shape[0]/2):,:]
                return self
            else:
                self.a = np.matmul(other, self.a)
                self.b = np.matmul(other, self.b)
                return self
        elif isinstance(other, complex):
            if hasattr(other, '__iter__'):
                if other.shape[1] == 2*self.a.shape[0]:
                    ab = other*np.concatenate([self.a, self.b])
                    return self
        else:
            #raise #TypeError('The multiplication operation is not possible for the type of variable used!')
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (int, float)):
            self.a = other*self.a
            self.b = other*self.b
            return self
        elif isinstance(other, np.ndarray):
            if other.shape[1] == 2*self.a.shape[0]:
                ab = np.matmul(other, np.concatenate([self.a, self.b]))
                self.a = ab[:int(ab.shape[0]/2),:]
                self.b = ab[int(ab.shape[0]/2):,:]
                return self
            else:
                self.a = np.matmul(other, self.a)
                self.b = np.matmul(other, self.b)
                return self
        elif isinstance(other, complex):
            if hasattr(other, '__iter__'):
                if other.shape[1] == 2*self.a.shape[0]:
                    ab = np.matmul(other, np.concatenate([self.a, self.b]))
                    self.a = ab[:int(ab.shape[0]/2),:]
                    self.b = ab[int(ab.shape[0]/2):,:]
                    return self
        else:
            #raise #TypeError('The multiplication operation is not possible for the type of variable used!')
            return NotImplemented
    __array_priority__ = 10000

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

    function bsc = beam(bsc, idx)
      % BEAM get beams from a beam array object
      %
      % BEAM(idx) idx can be a linear index or a logical array.
      bsc.a = bsc.a(:, idx);
      bsc.b = bsc.b(:, idx);
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