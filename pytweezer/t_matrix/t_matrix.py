import numpy as np
from pytweezer.shapes import *
from pytweezer.utils import combined_index

class TMatrix:
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
        self.k_medium = self.parser_k_medium(k_medium, 2.0*np.pi)    
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
        

    def parser_wavenumber(self, p, default):
        k_m = self.parser_k_medium(p, default=None)
        k_p = self.parser_k_particle(p, default=None)
        if not k_m and not k_p:
            k_m = default
        if p['index_r']:
            if not k_m and not k_p:
               km = k_p / p['index_r']
            elif k_m and not k_p:
                k_p = k_m * p['index_r']
            else:
                raise ValueError('index_relative specified but both indices already known');
        elif not k_p:
            raise ValueError('Unable to determine particle wavenumber from inputs');
        return k_m, k_p


    def parser_k_medium(self, parameters, default=None):
        if parameters.get('k_m'):
            return parameters.get('k_m')
        elif parameters.get('index_m'):
            if not parameters.get('lambda_0'):            
                raise ValueError('wavelength0 must be specified to use index_medium')
            return parameters.get('index_m')*2.0*np.pi/parameters.get('lambda_0')
        else:
            return default

    def parser_k_particle(self, parameters, default=None):
        if parameters.get('k_p'):
            return parameters['k_p']
        elif parameters.get('lambda_p'):
            return 2.0*np.pi/parameters['lambda_p']
        elif parameters.get('index_p'):
            if not parameters.get('lambda_0'):
                raise ValueError('wavelength0 must be specified to use index_particle')
            return parameters['index_p']*2*np.pi/parameters['lambda_0']
        else:
            return default

    def cylinder_preferred_method(self, parameters, k, tolerance=0.01):
        ebcm1 = {}
        ebcm1['x'] = np.array([73, 169, 198, 228, 261, 391, 586, 718, 718, 657, 523, 457, 262, 73])
        ebcm1['y'] = np.array([409, 406, 418, 423, 397, 412, 400, 375, 223, 193, 195, 165, 204, 390])
        ebcm1['x'] = (ebcm1['x'] - 73) * 2.0 / (718 - 73)
        ebcm1['y'] = -(ebcm1['y'] - 438) * 6.0 / (438 - 9)
        pm1 = {}
        pm1['x'] = np.array([297, 355, 394, 718, 718, 591, 525, 391, 361, 297])
        pm1['y'] = np.array([943, 933, 946, 894, 868, 846, 874, 864, 913, 913])
        pm1['x'] = (pm1['x'] - 73) * 2.0 / (718 - 73)
        pm1['y'] = -(pm1['y'] - 985) * 6.0 / (985 - 555)
        ebcm10 = {}
        ebcm10['x'] = np.array([73, 193, 718, 718, 525, 328, 229, 160, 73])
        ebcm10['y'] = np.array([430, 426, 381, 37, 94, 177, 214, 274, 375])
        ebcm10['x'] = (ebcm10['x'] - 73) * 2.0 / (718 - 73)
        ebcm10['y'] = -(ebcm10['y'] - 438) * 6.0 / (438 - 9)
        pm10 = {}
        pm10['x'] = np.array([130, 160, 328, 397, 462, 522, 589, 718, 718, 654, 589, 522, 328, 265, 130])
        pm10['y'] = np.array([961, 970, 967, 951, 946, 946, 925, 912, 753, 784, 798, 798, 865, 874, 948])
        pm10['x'] = (pm10['x'] - 73) * 2.0 / (718 - 73)
        pm10['y'] = -(pm10['y'] - 985) * 6.0 / (985 - 555)
        k = k * 1.064 / 2.0 / np.pi
        diameter = 2.0 * parameters[0] * k
        len = parameters[1] * k
        if tolerance >= 0.1:
            if inpolygon(diameter, len, pm10.x, pm10.y):
                method = 'pm'
            elif inpolygon(diameter, len, ebcm10.x, ebcm10.y):
                method = 'ebcm'
            else:
                method = 'other'
        elif tolerance >= 0.01:
            if inpolygon(diameter, len, pm1.x, pm1.y):
                method = 'pm'
            elif inpolygon(diameter, len, ebcm1.x, ebcm1.y):
                method = 'ebcm'
            else:
                method = 'other'
        else:
            method = 'other'

    def set_type(self, new_type):
        if not new_type in ['total', 'internal', 'scattered']:
            raise ValueError('Invalid T-Matrix type')
        if self.type:
            if self.type != new_type:
                if self.type == 'scattered' and new_type == 'total':
                    self.T = 2.0*self.T + np.eye(self.T.shape[0]) 
                elif self.type == 'total' and new_type == 'scattered':
                    self.T = 0.5*(self.T - np.eye(self.T.shape[0]))
                else:
                    raise ValueError('No known conversion')
        self.type = new_type

    def get_n_max(self):
        nmax1, _ = combined_index(self.T.shape[0]/2)
        nmax2, _ = combined_index(self.T.shape[1]/2)
        n_max = np.array([nmax1, nmax2]).astype(int)
        return n_max    
    
    def set_n_max(self, n_max, tolerance=1e-6, powerloss='warn'):
        if n_max.size == 2:
            n_max1, n_max2 = n_max
        else:
            n_max1 = n_max[0]
            n_max2 = n_max[0]
        if not all(np.array(n_max1, n_max2) == self.get_n_max()):
            total_orders1 = combined_index(n_max1, n_max1)[0]
            total_orders2 = combined_index(n_max2, n_max2)[0]
            midpoint1, midpoint2 = self.T.shape[0]/2
            old_type = tmatrix.type
            if total_orders1 > midpoint1 or total_orders2 > midpoint2 and old_type=='total':
                self.set_type('scattered')
            A11 = self.T[:midpoint1, :midpoint2]
            A12 = self.T[:midpoint1, midpoint2:]
            A21 = self.T[midpoint:, :midpoint2]
            A22 = self.T[midpoint:end,  midpoint:]
            if total_orders1 > midpoint1:
                row_index, col_index, a = find(A11)
                A11 = sparse(row_index,col_index,a,total_orders1,midpoint2)
                row_index, col_index, a = find(A12)
                A12 = sparse(row_index,col_index,a,total_orders1,midpoint2)
                row_index, col_index, a = find(A21)
                A21 = sparse(row_index,col_index,a,total_orders1,midpoint2)
                [row_index,col_index,a] = find(A22);
                A22 = sparse(row_index,col_index,a,total_orders1,midpoint2);
            elif total_orders1 < midpoint1:
                A11 = A11[1:total_orders1, :]
                A12 = A12[1:total_orders1, :]
                A21 = A21[1:total_orders1, :]
                A22 = A22[1:total_orders1, :]
            if total_orders2 > midpoint2:
                [row_index,col_index,a] = find(A11);
                A11 = sparse(row_index,col_index,a,total_orders1,total_orders2);
                [row_index,col_index,a] = find(A12);
                A12 = sparse(row_index,col_index,a,total_orders1,total_orders2);
                [row_index,col_index,a] = find(A21);
                A21 = sparse(row_index,col_index,a,total_orders1,total_orders2);
                [row_index,col_index,a] = find(A22);
                A22 = sparse(row_index,col_index,a,total_orders1,total_orders2);
            elif total_orders2 < midpoint2:
                A11 = A11[:, 1:total_orders2]
                A12 = A12[:, 1:total_orders2]
                A21 = A21[:, 1:total_orders2]
                A22 = A22[:, 1:total_orders2]
            if total_orders1 < midpoint1 or total_orders2 < midpoint2:
                magA = np.power(np.abs(self.T), 2).sum()
            tmatrix.T = np.array([[A11, A12], [A21, A22]])
            if total_orders1 < midpoint1 or total_orders2 < midpoint2:
                magA = np.power(np.abs(self.T), 2).sum()
                apparent_error = abs( magA - magB )/magA
                if apparent_error > tolerance:
                    if power_loss == 'warn':
                        warnings.warn('ott:Tmatrix:setNmax:truncation')
                    elif power_loss == 'error':
                        raise ValueError('ott:Tmatrix:setNmax:truncation')
                    elif power_loss=='ignore':
                        pass
                    else:
                        raise ValueError('Unrecognized option for powerloss')
            self.set_type(old_type)
        


'''
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
'''