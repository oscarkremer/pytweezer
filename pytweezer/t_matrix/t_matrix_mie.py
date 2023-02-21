
import ott.utils.sbesselj
import ott.utils.sbesselh1
from pytweezer.utils import combined_index, sbesselj, sbesselh

class TMatrixMie(TMatrix):
    def __init__(self):




'''
        p = ott.TmatrixMie.input_parser(varargin{:});
      tmatrix.radius = radiu
      [tmatrix.k_medium, tmatrix.k_particle] = ...
          tmatrix.parser_wavenumber(p, 2.0*pi);
      tmatrix.mu_relative = p.Results.mu_relative;
      if numel(tmatrix.radius) ~= numel(tmatrix.k_particle)
        error('radius and k_particle must be the same length');
      end
      if numel(tmatrix.radius) >= 100
        warning('May not work well for particles with >= 100 layers');
      end
      if isempty(p.Results.Nmax)
        if p.Results.internal == false
          Nmax = ott.utils.ka2nmax(tmatrix.k_medium*radius(end));
        else
          Nmax = ott.utils.ka2nmax(tmatrix.k_particle(end)*radius(end));
        end
      else
        Nmax = p.Results.Nmax;
      end
      if p.Results.internal == false
        tmatrix.type = 'scattered';
      else
        tmatrix.type = 'internal';
      end
      if numel(tmatrix.radius) == 1
        tmatrix.data = tmatrix.tmatrix_mie(Nmax, p.Results.internal);
      else
        if p.Results.shrink && isempty(p.Results.Nmax)
          oldNmax = Nmax;
          Nmax = max(100, Nmax);
        end
       tmatrix.data = tmatrix.tmatrix_mie_layered(Nmax, p.Results.internal);
        if p.Results.shrink
          tmatrix.Nmax = oldNmax;
        end

      end
    end
  end
end
        
        pass

    def tmatrix_mie(self, tmatrix, Nmax, internal):
        n = [1:Nmax]
        m = tmatrix.k_particle/tmatrix.k_medium
        mu = tmatrix.mu_relative
        r0 = tmatrix.k_medium * tmatrix.radius
        r1 = tmatrix.k_particle * tmatrix.radius

        indexing = combined_index(1:Nmax^2+2*Nmax)
        j0 = (sbesselj(n,r0)).';
        j1 = (sbesselj(n,r1)).';
        h0 = (sbesselh1(n,r0)).';
        j0d = (sbesselj(n-1,r0) - n.*sbesselj(n,r0)/r0).';
        j1d = (sbesselj(n-1,r1) - n.*sbesselj(n,r1)/r1).';
        h0d = (sbesselh1(n-1,r0) - n.*sbesselh1(n,r0)/r0).';

        if internal == false:
            b = -( mu*j1d.*j0 - m*j0d.*j1 ) ./ ( mu*j1d.*h0 - m*h0d.*j1 );
            a = -( mu*j0d.*j1 - m*j1d.*j0 ) ./ ( mu*h0d.*j1 - m*j1d.*h0 );
            T=sparse([1:2*(Nmax^2+2*Nmax)],[1:2*(Nmax^2+2*Nmax)], ...
                [a(indexing);b(indexing)]);
        else:
            d = ( h0d.*j0 - j0d.*h0 ) ./ ( m*j1.*h0d - j1d.*h0 );
            c = ( j0d.*h0 - h0d.*j0 ) ./ ( m*j1d.*h0 - h0d.*j1 );
            T=sparse([1:2*(Nmax^2+2*Nmax)],[1:2*(Nmax^2+2*Nmax)], ...
                [c(indexing);d(indexing)]);
        return T 

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