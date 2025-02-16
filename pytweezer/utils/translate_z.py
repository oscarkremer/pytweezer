import warnings
warnings.filterwarnings('ignore')
import numpy as np
from scipy.sparse import csr_matrix
from .sbesselj import sbesselj
from .sbesselh import sbesselh
from .combined_index import combined_index

def translate_z(nmax: int, z, function_type='sbesselj', method='gumerov'):
    #TODO docstirng
    if hasattr(z, '__iter__'):
        A, B, C = [], [], []
        for element in z:
            A_i, B_i, C_i = translate_z(nmax, element, function_type=function_type, method=method)
            A.append(A_i)
            B.append(B_i)
            C.append(C_i)
        return A, B, C

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

    C = C[:nmax+1, :nmax+1, :min(nmax1, nmax2)+1]
    return A, B, C

def translate_z_videen(nmax1, nmax2, nmax, z, function_type):
    N = 3*nmax+5
    N3 = min(nmax1, nmax2) + 1

    C = np.zeros((N,N,N3))+0j
    k = np.arange(0, N)
    if function_type not in ('sbesselj', 'sbesselh1', 'sbesselh2'):
        raise ValueError('Unknown value for function_type parameter, allowed values are: \
            (\'sbesselj\', \'sbesselh1\', \'sbesselh2\')')
    elif function_type == 'sbesselj':
        if z < 0:
            jn, _ = sbesselj(k, 2*np.pi*abs(z))
            C[:, 0, 0] = np.sqrt(2*k+1)*jn*(-1)**(k)
        else:
            jn, _ = sbesselj(k, 2*np.pi*z)
            C[:, 0, 0] = np.sqrt(2*k+1)*jn
    elif function_type == 'sbesselh1':
        if z < 0:
            hn, _ = sbesselh(k, 2*np.pi*abs(z), htype='1')
            C[:, 0, 0] = np.sqrt(2*k+1)*hn*(-1)**(k)
        else:
            hn, _ = sbesselh(k, 2*np.pi*z, htype='1')
            C[:, 0, 0] = np.sqrt(2*k+1)*hn
    elif function_type == 'sbesselh2':
        if z < 0:
            hn, _ = sbesselh(k, 2*np.pi*abs(z), htype='2')
            C[:, 0, 0] = np.sqrt(2*k+1)*hn*(-1)**(k)
        else:
            hn, _ = sbesselh(k, 2*np.pi*z, htype='2')
            C[:, 0, 0] = np.sqrt(2*k+1)*hn
    kk = np.arange(1, N-1)
    kk = kk.T.flatten()

    Cm = C[kk[0]-1:kk[-1], 0,0]
    Cm = Cm.T.flatten()
    Cp = C[kk[0]+1:kk[-1]+2, 0, 0]
    Cp = Cp.T.flatten()
    C[0, 1, 0] = -C[1, 0, 0]
    C[kk[0]:kk[-1]+1, 1, 0] = np.sqrt(3/(2*kk+1))*(kk*np.sqrt(1/(2*kk-1))*Cm-(kk+1)*np.sqrt(1/(2*kk+3))*Cp)

#    % Now do the rest, up to n=N-1
#    % Videen (40), with n(Videen) = n-1, n' = k
#    % Note that only the k=0 term is needed for n=N-1
    for i in range(2, N-1):
        kk = np.arange(1, N-i)
        kk = kk.T.flatten()
        Cm = C[kk[0]-1:kk[-1], i-1, 0]
        Cm = Cm.T.flatten()
        Cp = C[kk[0]+1:kk[-1]+2, i-1, 0]
        Cp = Cp.T.flatten()
        C0 = C[kk[0]:kk[-1]+1, i-2, 0].T.flatten()
        C[0, i, 0] = (-1)**i*C[i, 0, 0]
        C[kk[0]:kk[-1]+1, i, 0] = np.sqrt((2*i+1)/(2*kk+1))/i*(kk*np.sqrt((2*i-1)/(2*kk-1))*Cm+(i-1)*np.sqrt((2*kk+1)/(2*i-3))*C0-(kk+1)*np.sqrt((2*i-1)/(2*kk+3))*Cp)
    n = N-1
    C[0, N-1, 0] = np.sqrt(2*n+1)/n*((n-1)*np.sqrt(1/(2*n-3))*C[0, n-2, 0]-np.sqrt((2*n-1)/3)*C[1, n-1, 0])
#    % OK, now m other than m=0
#    % Only need to do positive m, since C(-m) = C(m)
#    % Videen (41)
    for m  in range(1, min(nmax1, nmax2)+1):
        nn = np.arange(m, nmax1+1)
        kk = np.arange(m, N-1)
        C0 = C[kk[0]:kk[-1]+1, nn[0]:nn[-1]+1, m-1]
        Cp = C[kk[0]+1:kk[-1]+2, nn[0]:nn[-1]+1,m-1]
        Cm = C[kk[0]-1:kk[-1], nn[0]:nn[-1]+1, m-1]
        factor1 = np.sqrt(1/((2*kk.reshape((-1,1))+1)*((nn.reshape((1,-1))-m+1)*(nn.reshape((1,-1))+m))))
        factor2 = np.sqrt(((kk-m+1)*(kk+m)*(2*kk+1))).reshape((-1,1))*C0
        factor3 = 2*np.pi*z*np.sqrt((((kk-m+2)*(kk-m+1)))/((2*kk+3))).reshape((-1, 1))*Cp
        factor4 = 2*np.pi*z*np.sqrt((((kk+m)*(kk+m-1)))/((2*kk-1))).reshape((-1, 1))*Cm
        C[kk[0]:kk[-1]+1, nn[0]:nn[-1]+1, m] = factor1*(factor2-factor3-factor4)
    return C

def calculate_AB(C, nmax1, nmax2, nmax, z, p):
    #    % OK, that's the scalar coefficients
    #% Time to find the vector coefficients - Videen (43) & (44)

    nn = np.arange(1, nmax1+1)
    kk = np.arange(1, nmax2+1).reshape((1,-1)).T
    matrixm = np.sqrt(kk*(kk+1))/np.sqrt(nn*(nn+1))

    central_iterator1 = np.arange(1, nmax1+1)*np.arange(2, nmax1+2)
    central_iterator2 = np.arange(1, nmax2+1)*np.arange(2, nmax2+2)
    ciy,cix = np.meshgrid(central_iterator1, central_iterator2);

    mmm = 0
    C0 = C[1:(nmax2+1), 1:(nmax1+1), mmm]
    Cp = C[2:(nmax2+2), 1:(nmax1+1), mmm]
    Cm = C[:nmax2, 1:(nmax1+1), mmm]
    t = matrixm*(C0 - 2*np.pi*np.abs(z)/(kk+1)*np.sqrt((kk-mmm+1)*(kk+mmm+1)/((2*kk+1)*(2*kk+3)))*Cp - 2*np.pi*np.abs(z)/kk*np.sqrt((kk-mmm)*(kk+mmm)/((2*kk+1)*(2*kk-1)))*Cm)
    toIndexy = ciy.T.flatten()
    toIndexx = cix.T.flatten()
    A = t.T.flatten()
    B = np.zeros(A.size)
    for mmm in range(1, min(nmax1, nmax2)+1):
        sz1 = np.arange(mmm, nmax2+1)
        sz2 = np.arange(mmm, nmax1+1)
        C0 = C[mmm:(nmax2+1), mmm:(nmax1+1), mmm]
        Cp = C[1+mmm:(nmax2+2), mmm:(nmax1+1), mmm]
        Cm = C[mmm-1:nmax2, mmm:(nmax1+1), mmm]
        tt = matrixm[sz1[0]-1:sz1[-1], sz2[0]-1:sz2[-1]]*(C0 - 2*np.pi*np.abs(z)/(kk[sz1-1]+1)*np.sqrt((kk[sz1-1]-mmm+1)*(kk[sz1-1]+mmm+1)/((2*kk[sz1-1]+1)*(2*kk[sz1-1]+3)))*Cp -2*np.pi*np.abs(z)/kk[sz1-1]*np.sqrt((kk[sz1-1]-mmm)*(kk[sz1-1]+mmm)/((2*kk[sz1-1]+1)*(2*kk[sz1-1]-1)))*Cm)
        ciys = ciy[mmm-1:, mmm-1:]
        cixs = cix[mmm-1:, mmm-1:]
        toIndexy = np.concatenate([toIndexy, (ciys.T.flatten()+mmm), (ciys.T.flatten()-mmm)])
        toIndexx = np.concatenate([toIndexx, (cixs.T.flatten()+mmm), (cixs.T.flatten()-mmm)])

        A = np.concatenate([A, tt.T.flatten(), tt.T.flatten()])
        tt = mmm/(kk[sz1-1]*(kk[sz1-1]+1))*matrixm[sz1[0]-1:sz1[-1], sz2[0]-1:sz2[-1]]* C0
        B = np.concatenate([B, tt.T.flatten(), -tt.T.flatten()])
    # Keep B real until the end, makes things run faster
    B = 1j*2*np.pi*np.abs(z)*B
    #% This is faster than A = A + sparse(...) and A(sub2ind(...)) = [...]
    if z < 0:
        n1, _ = combined_index(toIndexy)
        n2, _ = combined_index(toIndexx)
        B = csr_matrix((B*(-1)**(n1-n2+1), (toIndexy-1, toIndexx-1)), shape=(nmax1*(nmax1+2),nmax2*(nmax2+2))).toarray()
        A = csr_matrix((A*(-1)**(n1-n2), (toIndexy-1, toIndexx-1)), shape=(nmax1*(nmax1+2),nmax2*(nmax2+2))).toarray()
    else:
        B = csr_matrix((B, (toIndexy-1, toIndexx-1)), shape=(nmax1*(nmax1+2),nmax2*(nmax2+2))).toarray()
        A = csr_matrix((A, (toIndexy-1, toIndexx-1)), shape=(nmax1*(nmax1+2),nmax2*(nmax2+2))).toarray()
    return A, B

#@jit(nopython=True)
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
        shn, dhn = sbesselh(nd, kr, htype='1')
        C_nd00 = np.sqrt(2*nd+1)*shn/2
    elif function_type == 'sbesselh2':
        shn, _ = sbesselh(nd, kr, htype='2')
        C_nd00 = np.sqrt(2*nd+1)*shn/2
    C_ndn0 = np.zeros((nd.size + 1, nd.size + 1))+0j
    C_ndn0[1:C_nd00.size+1,1] = C_nd00
    C_ndn0[1,1:C_nd00.size+1] = (-1)**(nd)*C_nd00
    for j in range(1, nmax+1):
        i = np.arange(j, fval-j+1)
        num = anm_l(i[0]-2,0)*C_ndn0[i+1, i[0]-1]-anm_l(i,0)*C_ndn0[i+2, i[0]] + anm_l(i-1,0)*C_ndn0[i, i[0]]
        C_ndn0[i+1, i[0]+1] = num/anm_l(i[0]-1,0)

        C_ndn0[i[0]+1, i+1] = (-1)**(j+i)*C_ndn0[i+1, i[0]+1]
    C = np.zeros((nmax2+2, nmax1+1, mmax+1))+0j
    C[:, :, 0] = C_ndn0[1:(nmax2+3), 1:(nmax1+2)]
    ANM = anm_l(np.arange(0, 2*nmax+2).reshape((1, -1)).T,np.arange(1, nmax+1))
    IANM = 1/ANM

    for m in range(1, mmax+1):
        nd = np.arange(m, fval-m+1)
        C_nd1m = (bnm_l(nd, -m)*C_ndn0[nd, m]-bnm_l(nd+1, m-1)*C_ndn0[nd+2, m])/bnm_l(m, -m)
        C_ndn1 = np.zeros(C_ndn0.shape)+0j
        C_ndn1[m+1:C_nd1m.size+m+1, m+1] = C_nd1m
        C_ndn1[m+1, m+1:C_nd1m.size+m+1] = (-1)**(nd+m)*C_nd1m
        for j in range(m+1, nmax+1):
            i= np.arange(j, fval-j+1)
            sub_expression = ANM[i[0]-2, m-1]*C_ndn1[i+1, i[0]-1]-ANM[i,m-1]*C_ndn1[i+2,i[0]]+ANM[i-1,m-1]*C_ndn1[i,i[0]]
            C_ndn1[i+1, i[0]+1]= sub_expression*IANM[i[0]-1, m-1]
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