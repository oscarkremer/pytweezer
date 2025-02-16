import numpy as np
from numba import njit
from numba.pycc import CC

cc = CC('_nb_anm_bnm_')

#@cc.export('_nb_legendre_row_', 'Tuple((f8[:],f8[:,:]))(i8, f8[:])')
@njit(nopython=True, cache=True)
def bnm_l(n, m):
    b_nm = (2*(m<0)-1)*np.sqrt((n-m-1)*(n-m)/(2*n-1)/(2*n+1))
    m_bigger_n = np.argwhere(np.abs(m) > n)
    if m_bigger_n.size:
        m_bigger_n = tuple(m_bigger_n.T)
        b_nm[m_bigger_n] = 0
    return b_nm
  
if __name__=='__main__':
    cc.compile()

