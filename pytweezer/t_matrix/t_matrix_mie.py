
import warnings
import numpy as np
from scipy.sparse import csr_matrix
from .t_matrix import TMatrix
from pytweezer.utils import combined_index, sbesselj, sbesselh, ka_nmax

class TMatrixMie(TMatrix):
    def __init__(self, radius, **kwargs):
        parameters = self.input_parser(**kwargs)
        self.radius = radius
        if isinstance(self.radius, np.ndarray):
            if self.radius.size == 1:
                self.radius = self.radius[0]
        self.k_m, self.k_p = self.parser_wavenumber(parameters, 2.0*np.pi)
        self.mu_relative = parameters['mu_relative']

        if isinstance(self.radius, np.ndarray):
            if self.radius.size != self.k_p.size:
                raise ValueError('radius and k_particle must be the same length')    
            if self.radius.size >= 100:
                warnings.warn('May not work well for particles with >= 100 layers')
        if not parameters['n_max']:
            if not parameters['internal']:
                if isinstance(self.radius, np.ndarray):
                    n_max = ka_nmax(self.k_m*self.radius[-1])
                else:
                    n_max = ka_nmax(self.k_m*self.radius)
            else:
                n_max = ka_nmax(self.k_p[-1]*self.radius[-1])
        else:
            n_max = parameters['n_max']
        if not parameters['internal']:
            self.type = 'scattered'
        else:
            self.type = 'internal'
        if isinstance(self.radius, (int, float)): 
            self.T = self.t_matrix_mie(n_max, parameters['internal'])
        else:
            if parameters['shrink'] and not parameters['n_max']:
                old_nmax = n_max
                n_max = max(100, n_max)
            self.T = self.tmatrix_mie_layered(n_max, parameters['internal'])
            if parameters['shrink']:
                self._n_max_ = old_nmax

    def input_parser(self, **kwargs):
        if not kwargs.get('n_max'):
            kwargs['n_max'] = None
        if not kwargs.get('lambda_0'):
            kwargs['lambda_0'] = None
        if not kwargs.get('internal'):
            kwargs['internal'] = False
        if not kwargs.get('shrink'):
            kwargs['shrink'] = True
        if not kwargs.get('index_r'):
            kwargs['index_r'] = None
        if not kwargs.get('mu_relative'):
            kwargs['mu_relative'] = 1.0
        if not kwargs.get('k_m'):
            kwargs['k_m'] = None
        if not kwargs.get('lambda_m'):
            kwargs['lambda_m'] = None
        if not kwargs.get('index_m'):
            kwargs['index_m'] = None
        if not kwargs.get('k_p'):
            kwargs['k_p'] = None
        if not kwargs.get('lambda_p'):
            kwargs['lambda_p'] = None
        if not kwargs.get('index_p'):
            kwargs['index_p'] = None
        return kwargs

    def t_matrix_mie(self, n_max, internal):
        n = np.arange(1, n_max+1, 1)
        m = self.k_p/self.k_m
        mu = self.mu_relative
        r0 = self.k_m * self.radius
        r1 = self.k_p * self.radius
        indexing, _ = combined_index(np.arange(1, n_max**2+2*n_max+1,1))
        indexing = indexing.astype(int)-1
        j0 = sbesselj(n, r0)[0].T
        j1 = sbesselj(n, r1)[0].T
        h0 = sbesselh(n, r0, htype='1')[0].T

        j0d = (sbesselj(n-1,r0)[0] - n*sbesselj(n,r0)[0]/r0).T
        j1d = (sbesselj(n-1,r1)[0] - n*sbesselj(n,r1)[0]/r1).T
        h0d = (sbesselh(n-1,r0, htype='1')[0] - n*sbesselh(n,r0, htype='1')[0]/r0).T

        if not internal:
            b = -( mu*j1d*j0 - m*j0d*j1)/(mu*j1d*h0 - m*h0d*j1)
            a = -( mu*j0d*j1 - m*j1d*j0)/(mu*h0d*j1 - m*j1d*h0)
            n_index = np.arange(1,2*(n_max**2+2*n_max)+1,1).astype(int)-1
            T = csr_matrix((np.concatenate([a[indexing], b[indexing]]).reshape((-1)), 
                (n_index, n_index)), shape=(n_index.max()+1,n_index.max()+1)).toarray()
        else:
            d = ( h0d*j0 - j0d*h0 )/( m*j1*h0d - j1d*h0 )
            c = ( j0d*h0 - h0d*j0 )/( m*j1d*h0 - h0d*j1 )
            T=sparse(np.arange(1,2*(Nmax**2+2*Nmax)+1,1),np.arange(1,2*(Nmax**2+2*Nmax)+1,1),
                [c(indexing),d(indexing)])
        return T 

'''
    def tmatrix_mie_layered(self, tmatrix, Nmax, internal)
      %TMATRIX_MIE code from tmatrix_mie_layered.m

      k_layer=[tmatrix.k_particle,tmatrix.k_medium]; % Medium on outside
      radius=[tmatrix.radius,tmatrix.radius(end)]; % medium same as final layer
      n_layer=k_layer/2/pi;

      n = 1:Nmax; %n for nmax

      lastElement=length(k_layer);
      if length(tmatrix.k_particle)>1
          [jN,jNd] = sbesselj(n,[k_layer.*radius, ...
              k_layer(2:end).*radius(1:end-1)]);
          [hN,hNd] = sbesselh1(n,[k_layer.*radius, ...
              k_layer(2:end).*radius(1:end-1)]);
      else
          [jN,jNd] = sbesselj(n,k_layer(:).*radius(:));
          [hN,hNd] = sbesselh1(n,k_layer(:).*radius(:));
      end
      jN=jN.';
      hN=hN.';
      jNd=jNd.';
      hNd=hNd.';

      d1_N=jNd./jN;
      d3_N=hNd./hN;
      r_N=jN./hN;

      ha_0=1/n_layer(1);
      hb_0=1;

      for ii=1:length(tmatrix.k_particle)
          ha_n=ha_0;
          hb_n=hb_0;

          m_0 = n_layer(ii);

          d1_0=d1_N(:,ii);
          d3_0=d3_N(:,ii);
          r_0=r_N(:,ii);

          if ii>1
              m_n = n_layer(ii-1);

              d1_n=d1_N(:,lastElement+ii-1);
              d3_n=d3_N(:,lastElement+ii-1);
              r_n=r_N(:,lastElement+ii-1);
          else
              m_n=1;

              d1_n=0;
              d3_n=0;
              r_n=0;
          end

          g0=m_0*ha_n-m_n*d1_n;
          g1=m_0*ha_n-m_n*d3_n;
          g0b=m_n*hb_n-m_0*d1_n;
          g1b=m_n*hb_n-m_0*d3_n;
          q_0=r_n./r_0;
          ha_0=(g1.*d1_0-q_0.*g0.*d3_0)./(g1-q_0.*g0);
          hb_0=(g1b.*d1_0-q_0.*g0b.*d3_0)./(g1b-q_0.*g0b);

      end
      m_0 = n_layer(lastElement-1);
      m_1 = n_layer(lastElement);

      m=m_0/m_1;

      d1_1=d1_N(:,lastElement);
      d3_1=d3_N(:,lastElement);
      r_1=r_N(:,lastElement);

      % %TM/TE coeffs...
      al_1=r_1.*(ha_0-m*d1_1)./(ha_0-m*d3_1);
      bl_1=r_1.*(m*hb_0-d1_1)./(m*hb_0-d3_1);
      b = -al_1;
      a = -bl_1;
      indexing = combined_index(1:Nmax^2+2*Nmax)
      if not internal:
        T=sparse(1:2*(Nmax^2+2*Nmax),1:2*(Nmax^2+2*Nmax),[a(indexing);b(indexing)]);
      else:
        r_0=(jN(:,lastElement)./jN(:,lastElement-1));
        d = r_0.*(d3_1 - d1_1 )  ./ (ha_0-m*d3_1);
        c = r_0.*(d3_1 - d1_1 )  ./ (m*hb_0 - d3_1);
        warning('ott:tmatrix_mie_layered:internalcoefficientwarning', ...
            ['The internal coefficients are for the outermost layer only...' ...
             ' the real ones are only defined for each layer.']);
        T=sparse(1:2*(Nmax^2+2*Nmax),1:2*(Nmax^2+2*Nmax),[c(indexing);d(indexing)]);
    return T
'''