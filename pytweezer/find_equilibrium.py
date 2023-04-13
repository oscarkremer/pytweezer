import numpy as np
import warnings

def find_equilibrium(z, fz):

    if not isinstance(z, np.ndarray) or not isinstance(fz, np.ndarray):
        raise TypeError('Inputs z and fz must be numpy arrays.')
    if z.size != fz.size:
        raise ValueError('Number of element in z and fz must be equal')
    zero_index = np.where(fz<0)[0][0]
    if not zero_index:
        warnings.warn('Ignoring fz<0 entries at start of vector')
        zero_index_1 = np.where(fz>0)[0][0]
        zero_index = np.where(fz[zero_index_1:])[0][0]+zero_index_1 - 1
    zmin = z.min()
    zmax = z.max()
    z = 2 * (z - zmin) / (zmax - zmin) - 1
    if zero_index.size:
        zrange = np.arange(max(zero_index-2,1),min(zero_index+2,z.size)+1)
        pz = np.polyfit(z[zrange], fz[zrange], 3)
        root_z = np.roots(pz)
        dpz = np.array([3*pz[0],2*pz[1],1*pz[2]])
        real_z = root_z[root_z.imag==0][0]
        roots_of_sign = np.polyval(dpz,real_z)
        zeq = real_z[roots_of_sign<0]
        try:
            eq = zeq[abs(zeq-z[zero_index])==min(abs(zeq-z[zero_index]))]
            eq = (eq + 1)/2*(zmax - zmin) + zmin
            return eq
        except Exception as e:
            return []
    else:
        return []