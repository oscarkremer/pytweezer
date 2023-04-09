import numpy as np
from .t_matrix import TMatrixMie
from pytweezer.utils import combined_index
from numpy import matlib as matlib
from copy import copy

def force_torque(ibeam, sbeam, position=np.array([[],[],[]]), 
        rotation=np.array([[],[],[]]), coherent=False):
    fx, fy, fz, tx, ty, tz, sx, sy, sz= 0, 0, 0, 0, 0, 0, 0, 0, 0
    assert position.shape[0] == 3 or position.size == 0, 'Position must either be a empty array or 3xN array'
    assert rotation.shape[0] == 3 or rotation.size == 0, 'Rotation must either be a empty array or 3x3N array'
    if isinstance(sbeam, TMatrixMie):
        T = copy(sbeam)
        n_positions = max(1, position.shape[1])
        n_rotations = max(1, rotation.shape[1]/3)
        if n_positions != 1 and n_rotations != 1 and n+positions != n_rotations:
            raise ValueError('OTT:forcetorque:nlocations Number of positions/rotations should be equal or 1')
        n_locations = max(n_positions, n_rotations)
        T.set_type('scattered')   
        n_beams = ibeam.get_n_beams()
        if coherent:
            n_beams = 1
        f = np.zeros((3*T.T.size, n_locations, n_beams))
        t = np.zeros((3*T.T.size, n_locations, n_beams))
        s = np.zeros((3*T.T.size, n_locations, n_beams))
        for i in range(n_locations):
            if position.size:
                if n_positions == 1:
                    aux_position = position
                else:
                    aux_position = position[:,i]
            else:
                aux_position = np.array([[],[],[]])
            if rotation.size:
                if n_rotations == 1:
                    aux_rotation = rotation
                else:
                    aux_rotation = rotation[:,i]
            else:
                aux_rotation = np.array([[],[],[]])
            s_beam, t_beam = ibeam.scatter(T, position=aux_position, rotation=aux_rotation)
            if coherent:
                sbeam = s_beam.merge_beams()
                tbeam = t_beam.merge_beams()
            fl, tl, sl = _force_torque_(t_beam, s_beam)
            f[:, i, :] = fl.reshape((3*T.size, 1, n_beams))
            t[:, i, :] = tl.reshape((3*T.size, 1, n_beams))
            z[:, i, :] = sl.reshape((3*T.size, 1, n_beams))
        f = np.squeeze(f)
        t = np.squeeze(t)
        s = np.squeeze(s)
        fx = f[1*np.arange(1, n_particles+1), :]
        fy = f[2*np.arange(1, n_particles+1), :]
        fz = f[3*np.arange(1, n_particles+1), :]
        tx = t[1*np.arange(1, n_particles+1), :]
        ty = t[2*np.arange(1, n_particles+1), :]
        tz = t[3*np.arange(1, n_particles+1), :]
        sx = s[1*np.arange(1, n_particles+1), :]
        sy = s[2*np.arange(1, n_particles+1), :]
        sz = s[3*np.arange(1, n_particles+1), :]
        return fx, fy, fz, tx, ty, tz, sx, sy, sz


def _force_torque_(i_beam, s_beam):
    if (i_beam.get_n_beams() != s_beam.get_n_beams()) and i_beam.get_n_beams() != 1 and s_beam.get_n_beams() != 1:
        raise ValueError('Beam objects must contain same number of beams or 1 beam');
    if i_beam.get_n_max() > s_beam.get_n_max():
        s_beam.set_n_max(i_beam.get_n_max())# = i_beam.get_n_max()
    elif i_beam.get_n_max() < s_beam.get_n_max():
        i_beam.set_n_max(s_beam.get_n_max())# = s_beam.get_n_max()
    s_beam = s_beam.total_field(i_beam)
    print(s_beam.a)
    a, b = i_beam.get_coefficients()
    p, q = s_beam.get_coefficients()
    n, m = i_beam.get_mode_indices()

    nmax = i_beam.get_n_max()
    b = 1j*b
    q = 1j*q
    addv = np.zeros((2*nmax+3,1))
    at = [a, matlib.repmat(addv, 1, a.shape[1])]
    bt = [b, matlib.repmat(addv, 1, b.shape[1])]
    pt = [p, matlib.repmat(addv, 1, p.shape[1])]
    qt = [q, matlib.repmat(addv, 1, q.shape[1])]

    ci = combined_index(n,m)
    np1 = 2*n+2
    cinp1 = ci+np1
    cinp1mp1 = ci+np1+1
    cinp1mm1 = ci+np1-1
    cimp1 = ci+1
    kimp = (m>n-1)
    anp1 = at[cinp1, :]
    bnp1 = bt[cinp1, :]
    pnp1 = pt[cinp1, :]
    qnp1 = qt[cinp1, :]
    anp1mp1 = at[cinp1mp1, :]
    bnp1mp1 = bt[cinp1mp1, :]
    pnp1mp1 = pt[cinp1mp1, :]
    qnp1mp1 = qt[cinp1mp1, :]
    anp1mm1 = at[cinp1mm1, :]
    bnp1mm1 = bt[cinp1mm1, :]
    pnp1mm1 = pt[cinp1mm1, :]
    qnp1mm1 = qt[cinp1mm1, :]
    amp1 = at[cimp1, :]
    bmp1 = bt[cimp1, :]
    pmp1 = pt[cimp1, :]
    qmp1 = qt[cimp1, :]
    amp1[kimp, :] = 0
    bmp1[kimp, :] = 0
    pmp1[kimp, :] = 0
    qmp1[kimp, :] = 0 
    a=a[ci, :]
    b=b[ci, :]
    p=p[ci, :]
    q=q[ci, :]

    Az = m/n/(n+1)*imag(-(a)*conj(b)+conj(q)*(p))
    bz_part1 = 1/(n+1)*np.sqrt(n*(n-m+1)*(n+m+1)*(n+2)/(2*n+3)/(2*n+1))
    bz_part2 = imag(anp1*conj(a)+bnp1*conj(b)-(pnp1)*conj(p)(qnp1)*conj(q))

    Bz = bz_part1*bz_part2
    fz = 2*(Az.sum()+Bz.sum())
    Axy=1j/n/(n+1)*np.sqrt((n-m)*(n+m+1))*(conj(pmp1)*q - conj(amp1)*b - conj(qmp1)*p + conj(bmp1)*a)
    bxy_part1 = 1j/(n+1)*np.sqrt(n*(n+2))/np.sqrt((2*n+1)*(2*n+3))
    bxy_part2 = sqrt((n+m+1)*(n+m+2))*( p*conj(pnp1mp1) + q* conj(qnp1mp1) -a*conj(anp1mp1) -b*conj(bnp1mp1)) 
    bxy_part3 = sqrt((n-m+1)*(n-m+2))*(pnp1mm1*conj(p) + qnp1mm1*conj(q) - anp1mm1*conj(a) - bnp1mm1*conj(b))
    Bxy=bxy_part1*(bxy_part2+bxy_part3)
    fxy=sum(Axy+Bxy)
    fx=real(fxy)
    fy=imag(fxy)
    return fx,fy,fz