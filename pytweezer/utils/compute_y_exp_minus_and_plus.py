
import numpy as np
from numba import njit
from numba.pycc import CC

cc = CC('compute_y_exp_minus_and_plus')

@njit(nopython=True, cache=True)
@cc.export('compute_y_exp_minus_and_plus_n', 'Tuple((c16[:,:], c16[:], c16[:]))(i8[:], i8[:,:], f8[:], f8[:,:])')
def compute_y_exp_minus_and_plus_n(m, mv, phiM, pnm):
    pnm = pnm[np.abs(m), :]
    neg_part = pnm[m<0,:]
    dim = neg_part.shape[0]
    pnm[:dim, :] = (np.power(-1, -mv[m<0,:]))*neg_part
    pnm[dim:,:] = pnm[m>=0,:]
    Y =  pnm*np.exp(1j*mv*phiM)
    expplus = np.exp(1j*phiM)
    expminus = np.exp(-1j*phiM)
    return Y, expplus, expminus


@njit(nopython=True, cache=True)
@cc.export('compute_y_exp_n', '(i8[:], i8[:,:], f8[:,:], f8[:,:])')
def compute_y_exp_n(m, mv, phiM, pnm):
    pnm = pnm[np.abs(m), :]
    neg_part = pnm[m<0,:]
    dim = neg_part.shape[0]
    pnm[:dim, :] = (np.power(-1, -mv[m<0,:]))*neg_part
    pnm[dim:,:] = pnm[m>=0,:]
    return pnm*np.exp(1j*mv*phiM)
    

@njit(nopython=True)
@cc.export('compute_y_theta_phi', '(c16[:,:], c16[:,:], c16[:], c16[:], i8, i8[:], i8[:], i8[:,:], f8[:])')
def compute_y_theta_phi(Y, Y2, expplus, expminus, n, m, mi, mv, theta):
    ymplus = np.append(Y[1:,:], np.zeros((1,theta.shape[0])),axis=0)
    ymminus= np.empty(Y.shape, dtype=Y.dtype)
    ymminus[0,:] = np.zeros((1, theta.shape[0]))
    ymminus[1:,:] = Y[:-1,:]
    Ytheta =  np.sqrt((n-mv+1)*(n+mv))/2*expplus*ymminus  - np.sqrt((n-mv)*(n+mv+1))/2*expminus*ymplus
    ymplus = Y2[2:, :]
    ymminus = Y2[:-2,:]
    Yphi = 1j/2* np.sqrt((2*n+1)/(2*n+3))*(np.sqrt((n+mv+1)*(n+mv+2))*expminus*ymplus+np.sqrt((n-mv+1)*(n-mv+2))*expplus*ymminus)
    Y = Y[n+mi,:]
    Yphi = Yphi[n+mi,:]
    Ytheta = Ytheta[n+mi,:]
    return Y, Ytheta, Yphi


if __name__=='__main__':
    cc.compile()
