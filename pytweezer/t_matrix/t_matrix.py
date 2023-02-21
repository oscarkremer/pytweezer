from pytweezer.shapes import *

class TMatrix:
    '''    % Class representing the T-matrix of a scattering particle or lens.
    % This class can either be instantiated directly or used as a base
    % class for defining custom T-matrix types.
    %
    % This class is the base class for all other T-matrix object, you
    % should inherit from this class when defining your own T-matrix
    % creation methods. This class doesn't inherit from ``double`` or ``single``,
    % instead the internal array type can be set at creation allowing the
    % use of different data types such as ``sparse`` or ``gpuArray``.
    %
    % This class is not a handle class, therefore, when using the class
    % methods you need to store the resulting T-matrix output, for example::
    %
    %   tmatrix = ott.Tmatrix();
    %   new_tmatrix = tmatrix.scattered();
    %
    % Properties
    %   - data        -- The T-matrix this class encapsulates
    %   - type (enum) -- Type of T-matrix (total or scattered)
    %
    % Methods
    %   - total()     -- Convert to a total-field T-matrix
    %   - scattered() -- Convert to a scattered-field T-matrix
    %   - real        -- Extract real part of T-matrix
    %   - imag        -- Extract imaginary part of T-matrix
    %
    % Static methods
    %   - simple()    -- Construct a simple particle T-matrix
    %
    % See also Tmatrix, simple, :class:`+ott.TmatrixMie`.

    % This file is part of the optical tweezers toolbox.
    % See LICENSE.md for information about using/distributing this file.
    properties (SetAccess=protected)
        data          % The matrix this class encapsulates
    end

    properties (Access=private)
        type_         % Type of T-matrix (actual value)
    end

    properties (Dependent)
        Nmax          % Current size of T-matrix
        type          % Type of T-matrix (total or scattered)
    end
    '''

    def __init__(self):
        pass

    def default_method(self, shape, parameters=np.array([]),     
        method_tol=None,
        k_medium=None,
        index_m=None,
        lambda_0=None
        ):
        
        '''      % Determine the appropriate method for a particular shape
        % Returns one of 'mie', smarties', 'dda', 'ebcm', or 'pm'.
        %
        % DEFAULTMETHOD(shape) determine the default method to use for
        % a ott.shapes.Shape obejct.
        %
        % DEFAULTMETHOD(name, parameters) determine the default method
        % for a shape described by its name and parameters.
        %
        % Supported shape names [parameters]:
        %   'sphere'          Spherical (or layered sphere) [ radius ]
        %   'cylinder'        z-axis aligned cylinder [ radius height ]
        %   'ellipsoid'       Ellipsoid [ a b c]
        %   'superellipsoid'  Superellipsoid [ a b c e n ]
        %   'cone-tipped-cylinder'      [ radius height cone_height ]
        %   'cube'            Cube [ width ]
        %   'axisym'          Axis-symetric particle [ rho(:) z(:) ]
        '''
        self.k_medium = self.parser_k_medium(k_medium, , 2.0*np.pi)    
        if isinstance(shape, Sphere):
            #TODO insert info about isSphere for ellipsoid and superellipsoid
            return 'mie'
        elif isinstance(shape, Ellipsoid): 
            # TODO: Insert condition for isEllipsoid for SuperEllipsoid
            return 'smarties'
        elif isinstance(shape, SuperEllipsoid):
            return 'pm'
        elif isinstance(shape, Cube) or isinstance(shape, RectangularPrism):
            return 'pm'
        elif isinstance(shape, Cylinder) or isinstance(shape, AxiSymLerp):
            if isinstance(shape, Cylinder):
                parameters = np.array([ shape.radius, shape.height])
            elif isinstance(shape, AxiSymLerp):
                parameters = np.array([shape.rho.max(),shape.z.max() - shape.z.min()])
            return self.cylinder_preferred_method(parameters, k_medium, p.Results.method_tol)
            if method == 'other':
                return 'dda'
        else:
            raise ValueError('ott:Tmatrix:simple:no_shape Unsupported particle shape')
        
    def parser_k_medium(self, k_medium=None, index_m=None, lambda_0=None, default=2*np.pi):
        if k_medium:
            return k_medium
        elif index_m:
            if not lambda_0:            
                raise ValueError('wavelength0 must be specified to use index_medium')
            return index_m*2.0*np.pi/lambda_0
        else:
            return default
        
    def simple(self, shape, parameters=np.array([]), method=None, method_tol=None):
        p = inputParser
        p.KeepUnmatched = true
        p.addOptional('parameters', [])
        p.addParameter('method', '')
        p.addParameter('method_tol', [])
        p.parse(varargin{:});
        unmatched = ott.utils.unmatchedArgs(p);
        
        % Parse remaining inputs for k_medium
        pk = inputParser;
        pk.KeepUnmatched = true;
        pk.addParameter('k_medium', []);
        pk.addParameter('wavelength_medium', []);
        pk.addParameter('index_medium', []);
        pk.addParameter('wavelength0', []);
        pk.parse(unmatched{:});
        k_medium = self.parser_k_medium(pk, 2.0*np.pi)

        if ischar(shape) and ~isempty(p.Results.parameters):
            shape = ott.shapes.Shape.simple(shape, p.Results.parameters);
        elif ~isa(shape, 'ott.shapes.Shape') || ~isempty(p.Results.parameters)
           raise TypeError('Must input either Shape object or string and parameters');
        
        method = p.Results.method;
        if isempty(method)
            method = ott.Tmatrix.defaultMethod(shape, ...
                'method_tol', p.Results.method_tol, ...
                'k_medium', k_medium);
        if method=='mie':
            tmatrix = ott.TmatrixMie.simple(shape, unmatched{:});
        elif method == 'smarties':
            tmatrix = ott.TmatrixSmarties.simple(shape, unmatched{:});
        elif method == 'pm':
            tmatrix = TmatrixPm.simple(shape, unmatched{:});
        elif method == 'dda':
            tmatrix = ott.TmatrixDda.simple(shape, unmatched{:});
        elif method == 'ebcm':
            tmatrix = ott.TmatrixEbcm.simple(shape, unmatched{:});
        else:
            raise ValueError('Unsupported method specified')
        return tmatrix

    def cylinder_preferred_method(self, parameters, k, tolerance=0.01):
        ebcm1 = {}
        ebcm1['x'] = [73, 169, 198, 228, 261, 391, 586, 718, 718, 657, 523, 457, 262, 73]
        ebcm1['y'] = [409, 406, 418, 423, 397, 412, 400, 375, 223, 193, 195, 165, 204, 390]
        ebcm1['x'] = (ebcm1['x'] - 73) * 2.0 / (718 - 73)
        ebcm1['y'] = -(ebcm1['y'] - 438) * 6.0 / (438 - 9)
        pm1 = {}
        pm1.x = [297, 355, 394, 718, 718, 591, 525, 391, 361, 297]
        pm1.y = [943, 933, 946, 894, 868, 846, 874, 864, 913, 913]
        pm1.x = (pm1.x - 73) * 2.0 / (718 - 73)
        pm1.y = -(pm1.y - 985) * 6.0 / (985 - 555)
        ebcm10 = {}
        ebcm10.x = [73, 193, 718, 718, 525, 328, 229, 160, 73]
        ebcm10.y = [430, 426, 381, 37, 94, 177, 214, 274, 375]
        ebcm10.x = (ebcm10.x - 73) * 2.0 / (718 - 73)
        ebcm10.y = -(ebcm10.y - 438) * 6.0 / (438 - 9)
        pm10 = {}
        pm10.x = [130, 160, 328, 397, 462, 522, 589, 718, 718, 654, 589, 522, 328, 265, 130]
        pm10.y = [961, 970, 967, 951, 946, 946, 925, 912, 753, 784, 798, 798, 865, 874, 948]
        pm10.x = (pm10.x - 73) * 2.0 / (718 - 73)
        pm10.y = -(pm10.y - 985) * 6.0 / (985 - 555)
        k = k * 1.064 / 2.0 / pi;
        diameter = 2.0 * parameters[0] * k
        len = parameters[1] * k

      if tolerance >= 0.1
        if inpolygon(diameter, len, pm10.x, pm10.y)
            method = 'pm';
        elseif inpolygon(diameter, len, ebcm10.x, ebcm10.y)
            method = 'ebcm';
        else
            method = 'other';
        end
      elseif tolerance >= 0.01
        if inpolygon(diameter, len, pm1.x, pm1.y)
          method = 'pm';
        elseif inpolygon(diameter, len, ebcm1.x, ebcm1.y)
          method = 'ebcm';
        else
          method = 'other';
        end
      else
        method = 'other';
      end
    end

    def parser_wavenumber(p, default):
        km = ott.Tmatrix.parser_k_medium(p, []);
        kp = ott.Tmatrix.parser_k_particle(p, []);

        if isempty(km) && isempty(kp)
            km = default;
        if ~isempty(p.Results.index_relative)
            if isempty(km) && ~isempty(kp)
                km = kp ./ p.Results.index_relative;
            elif ~isempty(km) && isempty(kp)
                kp = km .* p.Results.index_relative;
            else
                raise ValueError('index_relative specified but both indices already known');
        elif isempty(kp)
            raise ValueError('Unable to determine particle wavenumber from inputs');
        return km, kp

   

    function k_particle = parser_k_particle(p, default)
      %PARSER_K_PARTICLE helper to get k_particle from a parser object

      if ~isempty(p.Results.k_particle)
        k_particle = p.Results.k_particle;
      elseif ~isempty(p.Results.wavelength_particle)
        k_particle = 2.0*pi/p.Results.wavelength_particle;
      elseif ~isempty(p.Results.index_particle)
        if isempty(p.Results.wavelength0)
          error('wavelength0 must be specified to use index_particle');
        end
        k_particle = p.Results.index_particle ...
            * 2.0*pi/p.Results.wavelength0;
      elseif nargin == 2
        k_particle = default;
      else
        error('Unable to determine k_particle from inputs');
      end
    end
  end

 methods

  
    function nmax = get.Nmax(tmatrix)
      %get.Nmax calculate Nmax from the current T-matrix data
      nmax1 = ott.utils.combined_index(size(tmatrix.data, 1)/2);
      nmax2 = ott.utils.combined_index(size(tmatrix.data, 2)/2);

      % Support non-square T-matrices
      nmax = [nmax1 nmax2];
    end

    function tmatrix = set.Nmax(tmatrix, nmax)
      %set.Nmax resizes the T-matrix
      tmatrix = tmatrix.set_Nmax(nmax);
    end

    function tmatrix = set_Nmax(tmatrix, nmax, varargin)
      % SET_NMAX resize the T-matrix, with additional options
      %
      % SET_NMAX(nmax) sets the T-matrix Nmax.  nmax should be a
      % scarar or 2 element vector with row/column Nmax.
      %
      % SET_NMAX(..., 'tolerance', tol) use tol as the warning error
      % level tolerance for resizing the beam.
      %
      % SET_NMAX(..., 'powerloss', mode) action to take if a power
      % loss is detected.  Can be 'ignore', 'warn' or 'error'.

      p = inputParser;
      p.addParameter('tolerance', 1.0e-6);
      p.addParameter('powerloss', 'warn');
      p.parse(varargin{:});

      % Convert the input to row/column sizes
      if length(nmax) == 2
        nmax1 = nmax(1);
        nmax2 = nmax(2);
      else
        nmax1 = nmax(1);
        nmax2 = nmax(1);
      end

      % Check if we need to do anything
      if all([nmax1, nmax2] == tmatrix.Nmax)
        return;
      end

      total_orders1 = ott.utils.combined_index(nmax1, nmax1);
      total_orders2 = ott.utils.combined_index(nmax2, nmax2);

      midpoint1 = size(tmatrix.data, 1)/2;
      midpoint2 = size(tmatrix.data, 2)/2;

      % The current resizing method only works for scattered fields
      old_type = tmatrix.type;
      if total_orders1 > midpoint1 || total_orders2 > midpoint2 ...
          && strcmpi(old_type, 'total')
        tmatrix.type = 'scattered';
      end

      % Split T-matrix into quadrants
      A11 = tmatrix.data(1:midpoint1, 1:midpoint2);
      A12 = tmatrix.data(1:midpoint1, (midpoint2+1):end);
      A21 = tmatrix.data((midpoint1+1):end, 1:midpoint2);
      A22 = tmatrix.data((midpoint1+1):end, (midpoint2+1):end);

      % Resize rows
      if total_orders1 > midpoint1

        [row_index,col_index,a] = find(A11);
        A11 = sparse(row_index,col_index,a,total_orders1,midpoint2);

        [row_index,col_index,a] = find(A12);
        A12 = sparse(row_index,col_index,a,total_orders1,midpoint2);

        [row_index,col_index,a] = find(A21);
        A21 = sparse(row_index,col_index,a,total_orders1,midpoint2);

        [row_index,col_index,a] = find(A22);
        A22 = sparse(row_index,col_index,a,total_orders1,midpoint2);

      elseif total_orders1 < midpoint1

        A11 = A11(1:total_orders1, :);
        A12 = A12(1:total_orders1, :);
        A21 = A21(1:total_orders1, :);
        A22 = A22(1:total_orders1, :);

      end

      % Resize cols
      if total_orders2 > midpoint2

        [row_index,col_index,a] = find(A11);
        A11 = sparse(row_index,col_index,a,total_orders1,total_orders2);

        [row_index,col_index,a] = find(A12);
        A12 = sparse(row_index,col_index,a,total_orders1,total_orders2);

        [row_index,col_index,a] = find(A21);
        A21 = sparse(row_index,col_index,a,total_orders1,total_orders2);

        [row_index,col_index,a] = find(A22);
        A22 = sparse(row_index,col_index,a,total_orders1,total_orders2);

      elseif total_orders2 < midpoint2

        A11 = A11(:, 1:total_orders2);
        A12 = A12(:, 1:total_orders2);
        A21 = A21(:, 1:total_orders2);
        A22 = A22(:, 1:total_orders2);

      end

      if total_orders1 < midpoint1 || total_orders2 < midpoint2
        magA = full(sum(sum(abs(tmatrix.data).^2)));
      end

      % Recombined T-matrix from quadrants
      tmatrix.data = [ A11 A12; A21 A22 ];

      if total_orders1 < midpoint1 || total_orders2 < midpoint2
        magB = full(sum(sum(abs(tmatrix.data).^2)));
        apparent_error = abs( magA - magB )/magA;

        if apparent_error > p.Results.tolerance
          if strcmpi(p.Results.powerloss, 'warn')
            ott.warning('ott:Tmatrix:setNmax:truncation', ...
                ['Apparent error of ' num2str(apparent_error)]);
          elseif strcmpi(p.Results.powerloss, 'error')
            error('ott:Tmatrix:setNmax:truncation', ...
                ['Apparent error of ' num2str(apparent_error)]);
          elseif strcmpi(p.Results.powerloss, 'ignore')
            % Nothing to do
          else
            error('Unrecognized option for powerloss');
          end
        end
      end

      % If we were originally total field, convert back
      tmatrix.type = old_type;
    end

    function type = get.type(tmatrix)
      % Get the T-matrix type
      type = tmatrix.type_;
    end

    function tmatrix = set.type(tmatrix, type)
      % Set the T-matrix type, converting if needed
      tmatrix = tmatrix.set_type(type);
    end

    function tmatrix = set_type(tmatrix, type, varargin)
      % SET_TYPE set T-matrix type with additional options

      p = inputParser;
      p.addParameter('convert', true);
      p.parse(varargin{:});

      % Check the type is valid
      if ~strcmpi(type, 'total') && ~strcmpi(type, 'internal') ...
          && ~strcmpi(type, 'scattered')
        error('Invalid T-matrix type');
      end

      % Do type conversions
      if p.Results.convert && ~isempty(tmatrix.type) ...
          && ~strcmpi(tmatrix.type, type)
        if strcmpi(tmatrix.type, 'scattered') && strcmpi(type, 'total')
          tmatrix.data = 2.0*tmatrix.data + speye(size(tmatrix.data));
        elseif strcmpi(tmatrix.type, 'total') && strcmpi(type, 'scattered')
          tmatrix.data = 0.5*(tmatrix.data - speye(size(tmatrix.data)));
        else
          error('No known conversion');
        end
      end

      % Set the type
      tmatrix.type_ = type;
    end

    function tmatrix = total(tmatrix)
      % TOTAL convert T-matrix to total
      tmatrix.type = 'total';
    end

    function tmatrix = scattered(tmatrix)
      % SCATTERED convert T-matrix to scattered
      tmatrix.type = 'scattered';
    end

    function sbeam = mtimes(tmatrix,ibeam)
      %MTIMES provide T-matrix multiplication overload
      if isa(ibeam, 'ott.Bsc')
        sbeam = ibeam.scatter(tmatrix);
      else
        % Provide default matrix multiplication
        sbeam = tmatrix;
        sbeam.data = sbeam.data * ibeam;
      end
    end
    
    function tmatrixs = uminus(tmatrix)
      %UMINUS unary minus for T-matrix
      tmatrixs = tmatrix;
      tmatrixs.data = -tmatrixs.data;
    end

    function tmatrixs = real(tmatrix)
      % Extract real part of T-matrix
      tmatrixs = tmatrix;
      tmatrixs.data = real(tmatrixs.data);
    end

    function tmatrixs = imag(tmatrix)
      % Extract imaginary part of T-matrix
      tmatrixs = tmatrix;
      tmatrixs.data = imag(tmatrixs.data);
    end

    function check = columncheck(tmatrix)
      % Check the power in each column (non-absorbing T-matrix check)
      
      if strcmpi(tmatrix.type, 'scattered')
        check = sum(abs(2.*tmatrix.data+eye(size(tmatrix.data))).^2, 1);
      else
        check = sum(abs(tmatrix.data).^2, 1);
      end
      
    end
  end
end
