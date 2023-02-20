from .beam import Beam
import numpy as np
from pytweezer.utils import spherical_harmonics

class PointMatch(Beam):

    def __init__(self):
        pass

    def bsc_far_field(self, nn, mm, e_field, theta, phi, 
        icm=np.array([[]]), zero_rejection_level=1e-8,
        invert_coefficient_matrix=False):

        if not icm.size:
            coefficient_matrix = np.zeros((e_field.size, 2*nn.size))
            for n in range(1, max(nn)+1):
                ci = np.where(nn == n)[0]
                _, dtY, dpY = spherical_harmonics(n, mm[ci-1], theta, phi)
                print(dtY.shape)
                print(dtY)
                coefficient_matrix[:,ci-1] = np.concatenate([dpY,-dtY]) * np.power(1j,(n+1)/np.sqrt(n*(n+1)))
                coefficient_matrix[:,ci+nn.size-1] = np.concatenate([dpY,-dtY])*np.power(1j, n/np.sqrt(n*(n+1)))
            if invert_coefficient_matrix:
                icm = pinv(coefficient_matrix)
            if p.Results.invert_coefficient_matrix:
                exp_coeffs = np.matmul(icm, e_field)
            else:
                exp_coeffs = np.matmul(np.linalg.inv(e_field), coefficient_matrix)
        else:
            assert icm.shape[1] == e_field.size, 'Number of cols in coefficient matrix must match length(e_field)'
            exp_coeffs = icm * e_field
        
        fa = exp_coeffs[:int(exp_coeffs.shape[0]/2),:]
        fb = expansion_coefficients[int(exp_coeffs.shape[0]/2):,:]
        if zero_rejection_level:
            pwr = np.power(np.abs(fa),2)+np.power(np.abs(fb),2)
            non_zero = pwr > zero_rejection_level*pwr.max()
            nn = nn[non_zero]
            mm = mm[non_zero]
            fa = fa[non_zero]
            fb = fb[non_zero]
        a, b = self.make_beam_vector(fa, fb, nn, mm) 
        return a, b
'''  
#s    def bsc_focalplane(self, nn, mm, e_field, kr, theta, phi, varargin):
           % point match beam coefficients around focal plane
        %
        % [a, b, cm] = bsc_focalplane(nn, mm, e_field, kr, theta, phi, ...)
        % a, b are the beam coefficients.  cm is the coefficient matrix.
        %
        % nn, mm are the mode indices to include in the coefficient matrix.
        %
        % e_field is a vector of E-field to points to match.
        % The format should be [ Ex(:); Ey(:); Ez(:) ]
        %
        % kr, theta, phi are the coordinates of the Efield values.
        % These should be vectors of the same as length(e_field)/3.
        %
        % Optional named arguments:
        %   zero_rejection_level   val   removes modes with less power than
        %       zero_rejection_level.  Default 1e-8.  Use [] to disable.
        %   coefficient_matrix     mat   Coefficient matrix to use.
        %       default [].
        p = inputParser;
        p.addParameter('inv_coefficient_matrix', []);
        p.addParameter('zero_rejection_level', 1e-8);
        p.addParameter('invert_coefficient_matrix', nargout == 3);
        p.parse(varargin{:});
        
        assert(length(e_field) == numel(e_field), ...
        'e_field must be N element vector');
        assert(numel(e_field)/3 == numel(kr), 'kr must be same size as e_field/3');
        assert(numel(e_field)/3 == numel(theta), 'theta must be same size as e_field/3');
        assert(numel(e_field)/3 == numel(phi), 'phi must be same size as e_field/3');

        #      % Generate coefficient matrix
        icm = p.Results.inv_coefficient_matrix;
        if isempty(icm)
            coefficient_matrix = zeros(length(e_field), length(nn));
            for n = 1:length(nn)

                # Find RgM, RgN as appropriate for each mode
                [M,N] = ott.utils.vswfcart(nn(n),mm(n),kr,theta,phi,3);
                if rem(nn(n),2) == 0:
                    % Even n
                    MN = [ M(:,1); M(:,2); M(:,3) ];
                else:
    #              % Odd n
                    MN = [ N(:,1); N(:,2); N(:,3) ];
                coefficient_matrix(:,n) = MN;
        #% Invert coefficient matrix for icm
            if nargout == 3 || p.Results.invert_coefficient_matrix
                icm = pinv(coefficient_matrix);
            if p.Results.invert_coefficient_matrix
                expansion_coefficients = icm * e_field;
            else
                expansion_coefficients = coefficient_matrix \ e_field;
        else
            assert(size(icm, 2) == length(e_field),
                'Number of rows in coefficient matrix must match length(e_field)');
            expansion_coefficients = icm * e_field
       
#      % Look for non-zero elements, only keep non-zeros
        if ~isempty(p.Results.zero_rejection_level)
            non_zero = abs(expansion_coefficients) ...
                > max(abs(expansion_coefficients)) * p.Results.zero_rejection_level;
            expansion_coefficients = expansion_coefficients(non_zero);
            nn = nn(non_zero)
            mm = mm(non_zero)

#      % Calculate beam vectors
      fa = zeros(size(nn));
      fb = zeros(size(nn));
      for n = 1:length(nn):
         if rem(nn(n),2) == 0:
            fa(n) = expansion_coefficients(n);
            fb(n) = expansion_coefficients(n) * sign(mm(n));
         else:
            fa(n) = expansion_coefficients(n) * sign(mm(n));
            fb(n) = expansion_coefficients(n);
      [a, b] = ott.Bsc.make_beam_vector(fa, fb, nn, mm);

    def BscPointmatch(self, *args, **kwargs):
              % Protected constructor for BscPointmatch object
          %
          % See also ott.BscPmGauss and ott.BscPmParaxial
        return Bsc()

    def cleanCoefficientMatrix(self):
        # Remove the coefficient matrix data
        self.inv_coefficient_matrix = []


    
    function [a, b, cm] = bsc_focalplane(nn, mm, e_field, kr, theta, phi, ...
        varargin)
      % point match beam coefficients around focal plane
      %
      % [a, b, cm] = bsc_focalplane(nn, mm, e_field, kr, theta, phi, ...)
      % a, b are the beam coefficients.  cm is the coefficient matrix.
      %
      % nn, mm are the mode indices to include in the coefficient matrix.
      %
      % e_field is a vector of E-field to points to match.
      % The format should be [ Ex(:); Ey(:); Ez(:) ]
      %
      % kr, theta, phi are the coordinates of the Efield values.
      % These should be vectors of the same as length(e_field)/3.
      %
      % Optional named arguments:
      %   zero_rejection_level   val   removes modes with less power than
      %       zero_rejection_level.  Default 1e-8.  Use [] to disable.
      %   coefficient_matrix     mat   Coefficient matrix to use.
      %       default [].

      p = inputParser;
      p.addParameter('inv_coefficient_matrix', []);
      p.addParameter('zero_rejection_level', 1e-8);
      p.addParameter('invert_coefficient_matrix', nargout == 3);
      p.parse(varargin{:});
      
      assert(length(e_field) == numel(e_field), ...
        'e_field must be N element vector');
      assert(numel(e_field)/3 == numel(kr), 'kr must be same size as e_field/3');
      assert(numel(e_field)/3 == numel(theta), 'theta must be same size as e_field/3');
      assert(numel(e_field)/3 == numel(phi), 'phi must be same size as e_field/3');

      % Generate coefficient matrix
      icm = p.Results.inv_coefficient_matrix;
      if isempty(icm)
        coefficient_matrix = zeros(length(e_field), length(nn));
        for n = 1:length(nn)

           % Find RgM, RgN as appropriate for each mode
           [M,N] = ott.utils.vswfcart(nn(n),mm(n),kr,theta,phi,3);
           if rem(nn(n),2) == 0
              % Even n
              MN = [ M(:,1); M(:,2); M(:,3) ];
           else
              % Odd n
              MN = [ N(:,1); N(:,2); N(:,3) ];
           end
           coefficient_matrix(:,n) = MN;

        end

        % Invert coefficient matrix for icm
        if nargout == 3 || p.Results.invert_coefficient_matrix
          icm = pinv(coefficient_matrix);
        end

        % Do point matching
        if p.Results.invert_coefficient_matrix
          expansion_coefficients = icm * e_field;
        else
          expansion_coefficients = coefficient_matrix \ e_field;
        end
      else
        assert(size(icm, 2) == length(e_field), ...
            'Number of rows in coefficient matrix must match length(e_field)');

        % Do point matching
        expansion_coefficients = icm * e_field;
      end

      % Look for non-zero elements, only keep non-zeros
      if ~isempty(p.Results.zero_rejection_level)
        non_zero = abs(expansion_coefficients) ...
            > max(abs(expansion_coefficients)) * p.Results.zero_rejection_level;
        expansion_coefficients = expansion_coefficients(non_zero);
        nn = nn(non_zero);
        mm = mm(non_zero);
      end

      % Calculate beam vectors
      fa = zeros(size(nn));
      fb = zeros(size(nn));
      for n = 1:length(nn)

         if rem(nn(n),2) == 0
            fa(n) = expansion_coefficients(n);
            fb(n) = expansion_coefficients(n) * sign(mm(n));
         else
            fa(n) = expansion_coefficients(n) * sign(mm(n));
            fb(n) = expansion_coefficients(n);
         end

      end

      % Make the beam vector and store the coefficients
      [a, b] = ott.Bsc.make_beam_vector(fa, fb, nn, mm);
    end
  end
        '''
