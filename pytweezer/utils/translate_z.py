import warnings
import numpy as np
from pytweezer.utils import sbesselj

def translate_z(nmax: int, z, function_type='sbesselj', method='gumerov'):
    #TODO docstirng
    #    if hasattr(z, '__iter__'):
    #    A= np.cell(numel(z),1)
    #    B= A
    #    for element in z:
    #        A{ii},B{ii} =translate_z(nmax, element, function_type=function_type, method=method)
    #    C=0
#        ott.warning('external')
    #    return A, B, C

    #% Calculate Nmax for each dimension
    if hasattr(nmax, '__iter__'):
        if len(nmax) > 2:
            warnings.warn('Variable nmax given with more than two parameters, \
                using only the two first elements: {} and {}'.format(nmax[0], nmax[1]))
        nmax1 = nmax[0]
        nmax2 = nmax[1]
        nmax = max(nmax[0], nmax[1])
    else:
        nmax1 = nmax
        nmax2 = nmax    
    if z==0:
        A = np.eye(nmax1**2+nmax1*2, nmax2**2+nmax2*2)
        B = np.zeros((nmax1**2+nmax1*2, nmax2**2+nmax2*2))
        C = A
        return A, B, C

    #% Calculate the scalar coefficients
    if method not in ('videen', 'gumerov'):
        raise ValueError('Unknown method. Method for translation must be in (\'videen\',\'gumerov\'))')
    elif method == 'videen':
        C = translate_z_videen(nmax1, nmax2, nmax, abs(z), function_type)
    elif method == 'gumerov':
        C = translate_z_gumerov(nmax1, nmax2, nmax, abs(z), function_type)
 
    A, B = calculate_AB(C, nmax1, nmax2, nmax, z, function_type)

 #   C = C[:nmax+1, 1:nmax+1, 1:min(nmax1, nmax2)+1]
    return A, B, C

'''    
def translate_z_videen(nmax1, nmax2, nmax, z, p):
    N = 3*nmax+5;
    N3 = min(nmax1, nmax2) + 1;

    C = zeros(N,N,N3);

    #% First calculate the scalar translation coeffs

    #% Starting values, for m=0 and k=any -> n=0
    #% Videen (38)
    k = 0:(N-1);

switch p.Results.type
  case 'sbesselj'       % regular to regular
    if z < 0
      C(:,1,1) = sqrt(2*k+1) .* sbesselj(k,2*pi*abs(z)) .* (-1).^(k);
    else
      C(:,1,1) = sqrt(2*k+1) .* sbesselj(k,2*pi*z);
    end

  case 'sbesselh1'      % outgoing to regular
    if z < 0
      C(:,1,1) = sqrt(2*k+1) .* sbesselh1(k,2*pi*abs(z)) .* (-1).^(k);
    else
      C(:,1,1) = sqrt(2*k+1) .* sbesselh1(k,2*pi*z);
    end

  case 'sbesselh2'      % incoming to regular
    if z < 0
      C(:,1,1) = sqrt(2*k+1) .* sbesselh2(k,2*pi*abs(z)) .* (-1).^(k);
    else
      C(:,1,1) = sqrt(2*k+1) .* sbesselh2(k,2*pi*z);
    end

  case 'sbesselh1farfield'    % outgoing to regular

    if 2*pi*abs(z) <= (N-1).^2
      ott.warning('Farfield limit may not be satisfied');
    end

    h1limit = (-1i).^k./(1i*2*pi*abs(z)) .* exp(1i*2*pi*abs(z));

    if z < 0
      C(:,1,1) = sqrt(2*k+1) .* h1limit .* (-1).^(k);
    else
      C(:,1,1) = sqrt(2*k+1) .* h1limit;
    end

  case 'sbesselh2farfield'    % incoming to regular

    if 2*pi*abs(z) <= (N-1).^2
      ott.warning('Farfield limit may not be satisfied');
    end

    h2limit = (1i).^k./(-1i*2*pi*abs(z)) .* exp(-1i*2*pi*abs(z));

    if z < 0
      C(:,1,1) = sqrt(2*k+1) .* h2limit .* (-1).^(k);
    else
      C(:,1,1) = sqrt(2*k+1) .* h2limit;
    end

  otherwise
    error('OTT:UTILS:translate_z:type_error', 'Unknown translation type');
end

% Do n=1 as a special case (Videen (40) with n=0,n'=k)
kk = 1:(N-2);
kk = kk(:);
Cm = C(kk,1,1);
Cm = Cm(:);
Cp = C(kk+2,1,1);
Cp = Cp(:);
C(1,2,1) = -C(2,1,1);
C(kk+1,2,1) = sqrt(3./(2*kk+1)) .* ...
    ( kk.*sqrt(1./(2*kk-1)) .* Cm - (kk+1).*sqrt(1./(2*kk+3)) .* Cp );

% Now do the rest, up to n=N-1
% Videen (40), with n(Videen) = n-1, n' = k
% Note that only the k=0 term is needed for n=N-1
for n = 2:(N-2)
    kk = 1:(N-n-1);
    kk = kk(:);
    Cm = C(kk,n,1);
    Cm = Cm(:);
    Cp = C(kk+2,n,1);
    Cp = Cp(:);
    C0 = C(kk+1,n-1,1);
    C0 = C0(:);
    C(1,n+1,1) = (-1)^n * C(n+1,1,1);
    C(kk+1,n+1,1) = sqrt((2*n+1)./(2*kk+1))/n .* ...
        ( kk.*sqrt((2*n-1)./(2*kk-1)) .* Cm ...
        + (n-1)*sqrt((2*kk+1)/(2*n-3)) .* C0 ...
        - (kk+1).*sqrt((2*n-1)./(2*kk+3)) .* Cp );
end
n = N-1;
C(1,N,1) = sqrt(2*n+1)/n * ...
    ( (n-1)*sqrt(1/(2*n-3)) * C(1,n-1,1) - sqrt((2*n-1)/3) * C(2,n,1) );

% OK, now m other than m=0
% Only need to do positive m, since C(-m) = C(m)
% Videen (41)
for m = 1:min(nmax1, nmax2)
  nn = m:nmax1;
  kk = (m:N-2).';
  C0 = C(kk+1,nn+1,m);
  Cp = C(kk+2,nn+1,m);
  Cm = C(kk,nn+1,m);
  C(kk+1,nn+1,m+1) = sqrt(1./((2*kk+1)*((nn-m+1).*(nn+m)))) .* ...
      ( sqrt(((kk-m+1).*(kk+m).*(2*kk+1))).*C0 ...
      -2*pi*z*sqrt((((kk-m+2).*(kk-m+1)))./((2*kk+3))).*Cp ...
      -2*pi*z*sqrt((((kk+m).*(kk+m-1)))./((2*kk-1))).*Cm );
end

end % translate_z_videen
'''
def calculate_AB(C, nmax1, nmax2, nmax, z, p):
    #    % OK, that's the scalar coefficients
    #% Time to find the vector coefficients - Videen (43) & (44)

    nn = np.arange(1, nmax1+1)
    kk = np.arange(1, nmax2+1).reshape((1,-1)).T
    matrixm = np.sqrt(kk*(kk+1))/np.sqrt(nn*(nn+1))

    central_iterator1 = np.arange(1, nmax1+1)*np.arange(2, nmax1+2)
    central_iterator2 = np.arange(1, nmax2+1)*np.arange(2, nmax2+2)
    ciy,cix = np.meshgrid(central_iterator1,central_iterator2);

    mmm = 0
    C0 = C[1:(nmax2+1), 1:(nmax1+1), mmm]
    Cp = C[2:(nmax2+2), 1:(nmax1+1), mmm]
    Cm = C[:nmax2, 1:(nmax1+1), mmm]
    t = matrixm*(C0 - 2*np.pi*np.abs(z)/(kk+1)*np.sqrt((kk-mmm+1)*(kk+mmm+1)/((2*kk+1)*(2*kk+3)))*Cp - 2*np.pi*np.abs(z)/kk*np.sqrt((kk-mmm)*(kk+mmm)/((2*kk+1)*(2*kk-1)))*Cm)
    toIndexy=(ciy(:));
    toIndexx=(cix(:));
    A=t(:);
    B=zeros(size(A));


'''
for mmm=1:min(nmax1, nmax2)

    sz1 = mmm:nmax2;
    sz2 = mmm:nmax1;

    C0 = C((1+mmm):(nmax2+1),(1+mmm):(nmax1+1),mmm+1);
    Cp = C((2+mmm):(nmax2+2),(1+mmm):(nmax1+1),mmm+1);
    Cm = C((mmm):nmax2,(1+mmm):(nmax1+1),mmm+1);

    tt = matrixm(sz1, sz2).*(C0 - 2*pi*abs(z)./(kk(sz1)+1) .* ...
        sqrt((kk(sz1)-mmm+1).*(kk(sz1)+mmm+1) ...
        ./((2*kk(sz1)+1).*(2*kk(sz1)+3))) .* Cp - ...
        2*pi*abs(z)./kk(sz1).*sqrt((kk(sz1)-mmm) ...
        .*(kk(sz1)+mmm)./((2*kk(sz1)+1).*(2*kk(sz1)-1))).*Cm);

    ciys=ciy(mmm:end,mmm:end);
    cixs=cix(mmm:end,mmm:end);

    toIndexy=[toIndexy;(ciys(:)+mmm);(ciys(:)-mmm)];
    toIndexx=[toIndexx;(cixs(:)+mmm);(cixs(:)-mmm)];
    A=[A;tt(:);tt(:)];

    tt = mmm./(kk(sz1).*(kk(sz1)+1)).*matrixm(sz1, sz2) .* C0;
    B=[B;tt(:);-tt(:)];

end

% Keep B real until the end, makes things run faster
B = 1i*2*pi*abs(z)*B;

% This is faster than A = A + sparse(...) and A(sub2ind(...)) = [...]
if z < 0
  [n1, ~] = ott.utils.combined_index(toIndexy);
  [n2, ~] = ott.utils.combined_index(toIndexx);
  B=sparse(toIndexy,toIndexx,B.*(-1).^(n1-n2+1),nmax1*(nmax1+2),nmax2*(nmax2+2));
  A=sparse(toIndexy,toIndexx,A.*(-1).^(n1-n2),nmax1*(nmax1+2),nmax2*(nmax2+2));
else
  B=sparse(toIndexy,toIndexx,B,nmax1*(nmax1+2),nmax2*(nmax2+2));
  A=sparse(toIndexy,toIndexx,A,nmax1*(nmax1+2),nmax2*(nmax2+2));
end

end % calculate_AB
'''

def translate_z_gumerov(nmax1, nmax2, nmax, r, function_type):
    mmax = min(nmax1, nmax2)
    m = 0
    fval = 2*nmax+1
    nd = np.arange(m, fval+1)
    kr=2*np.pi*r
    if function_type not in ('sbesselj', 'sbesselh1', 'sbesselh2'):
        raise ValueError('Unknown value for function_type parameter, allowed values are: \
            (\'sbesselj\', \'sbesselh1\', \'sbesselh2\')')
    elif function_type == 'sbesselj':
        jn, _ = sbesselj(nd,kr) 
        C_nd00 = np.sqrt(2*nd+1)*jn
    elif function_type == 'sbesselh1':
        C_nd00=[np.sqrt(2*nd+1)*sbesselh(nd,kr, htype='1')]/2
    elif function_type == 'sbesselh2':
        C_nd00=[np.sqrt(2*nd+1)*sbesselh(nd, kr, htype='2')]/2
    C_ndn0 = np.zeros((nd.size + 1, nd.size + 1))
    C_ndn0[1:C_nd00.size+1,1] = C_nd00
    C_ndn0[1,1:C_nd00.size+1] = (-1)**(nd)*C_nd00
    for j in range(1, nmax+1):
        i = np.arange(j, fval-j+1)
        num = anm_l(i[0]-2,0)*C_ndn0[i+1, i[0]-1]-anm_l(i,0)*C_ndn0[i+2, i[0]] + anm_l(i-1,0)*C_ndn0[i, i[0]]
        C_ndn0[i+1, i[0]+1] = num/anm_l(i[0]-1,0)

        C_ndn0[i[0]+1, i+1] = (-1)**(j+i)*C_ndn0[i+1, i[0]+1]
    C = np.zeros((nmax2+2, nmax1+1, mmax+1))
    C[:, :, 0] = C_ndn0[1:(nmax2+3), 1:(nmax1+2)]
    ANM = anm_l(np.arange(0, 2*nmax+2).reshape((1, -1)).T,np.arange(1, nmax+1))
    IANM = 1/ANM
    for m in range(1, mmax+1):
        nd = np.arange(m, fval-m+1)
        C_nd1m = (bnm_l(nd, -m)*C_ndn0[nd, m]-bnm_l(nd+1, m-1)*C_ndn0[nd+2, m])/bnm_l(m, -m)
        C_ndn1 = np.zeros(C_ndn0.shape)
        C_ndn1[m+1:C_nd1m.size+m+1, m+1] = C_nd1m

        C_ndn1[m+1, m+1:C_nd1m.size+m+1] = (-1)**(nd+m)*C_nd1m

        for j in range(m+1, nmax+1):
            i= np.arange(j, fval-j+1)
            sub_expression = ANM[i[0]-2, m-1]*C_ndn1[i+1, i[0]-1]-ANM[i,m-1]*C_ndn1[i+2,i[0]]+ANM[i-1,m-1]*C_ndn1[i,i[0]]
            C_ndn1[i+1, i[0]+1]= sub_expression*IANM[i[0], m]
            C_ndn1[i[0]+1, i+1] = (-1)**(j+i)*C_ndn1[i+1,i[0]+1]
        C_ndn0=C_ndn1
        C[:,:,m] = C_ndn0[1:nmax2+3, 1:nmax1+2]
    return C



def anm_l(n, m):
    fn = 1/(2*n+1)/(2*n+3)
    a_nm = np.sqrt((n+np.abs(m)+1)*(n-np.abs(m)+1)*fn)
    smaller_than_zero = np.argwhere(n < 0)
    if smaller_than_zero.size:
        smaller_than_zero = tuple(smaller_than_zero.T)
        a_nm[np.argwhere(n < 0)] = 0
    m_bigger_n = np.argwhere(np.abs(m) > n)
    if m_bigger_n.size:
        m_bigger_n = tuple(m_bigger_n.T)
        a_nm[m_bigger_n] = 0
    return a_nm


def bnm_l(n, m):
    b_nm = (2*(m<0)-1)*np.sqrt((n-m-1)*(n-m)/(2*n-1)/(2*n+1))
    m_bigger_n = np.argwhere(np.abs(m) > n)
    if m_bigger_n.size:
        m_bigger_n = tuple(m_bigger_n.T)
        b_nm[m_bigger_n] = 0
    return b_nm