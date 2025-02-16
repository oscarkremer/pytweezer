import numpy as np
from .t_matrix import TMatrixMie
from pytweezer.utils import combined_index
from numpy import matlib as matlib
from copy import copy
import time


def force_torque(ibeam, sbeam, position=np.array([[],[],[]]), 
        rotation=np.array([[],[],[]]), coherent=False):
    fx, fy, fz, tx, ty, tz, sx, sy, sz= 0, 0, 0, 0, 0, 0, 0, 0, 0
    assert position.shape[0] == 3 or position.size == 0, 'Position must either be a empty array or 3xN array'
    assert rotation.shape[0] == 3 or rotation.size == 0, 'Rotation must either be a empty array or 3x3N array'
    if isinstance(sbeam, TMatrixMie):
        T = copy(sbeam)
        n_particles = 1
        n_positions = max(1, position.shape[1])
        n_rotations = max(1, rotation.shape[1]/3)
        if n_positions != 1 and n_rotations != 1 and n+positions != n_rotations:
            raise ValueError('OTT:forcetorque:nlocations Number of positions/rotations should be equal or 1')
        n_locations = max(n_positions, n_rotations)
        T.set_type('scattered')   
        n_beams = ibeam.get_n_beams()
        if coherent:
            n_beams = 1
        f = np.zeros((3, n_locations, n_beams))
        t = np.zeros((3, n_locations, n_beams))
        s = np.zeros((3, n_locations, n_beams))
        total_time = 0
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
            start = time.time()
            s_beam, t_beam = ibeam.scatter(T, position=aux_position, rotation=aux_rotation)
            end = time.time()
            part_time = end-start
            total_time+=part_time
            if coherent:
                sbeam = s_beam.merge_beams()
                tbeam = t_beam.merge_beams()
            fl = _force_torque_(t_beam, s_beam)
            f[:, i, :] = fl.reshape((3, n_beams))
            #t[:, i, :] = tl.reshape((3*T, 1, n_beams))
            #z[:, i, :] = sl.reshape((3, 1, n_beams))
        f = np.squeeze(f)
        #t = np.squeeze(t)
        #s = np.squeeze(s)
        fx = f[1*np.arange(1, n_particles+1)-1, :]
        fy = f[2*np.arange(1, n_particles+1)-1, :]
        fz = f[3*np.arange(1, n_particles+1)-1, :]
       # tx = t[1*np.arange(1, n_particles+1), :]
       # ty = t[2*np.arange(1, n_particles+1), :]
       # tz = t[3*np.arange(1, n_particles+1), :]
       # sx = s[1*np.arange(1, n_particles+1), :]
       # sy = s[2*np.arange(1, n_particles+1), :]
       # sz = s[3*np.arange(1, n_particles+1), :]
        return np.array([fx,fy,fz])# tx, ty, tz, sx, sy, sz


def _force_torque_(i_beam, s_beam):
    if (i_beam.get_n_beams() != s_beam.get_n_beams()) and i_beam.get_n_beams() != 1 and s_beam.get_n_beams() != 1:
        raise ValueError('Beam objects must contain same number of beams or 1 beam');
    if i_beam.get_n_max() > s_beam.get_n_max():
        s_beam.set_n_max(i_beam.get_n_max())# = i_beam.get_n_max()
    elif i_beam.get_n_max() < s_beam.get_n_max():
        i_beam.set_n_max(s_beam.get_n_max())# = s_beam.get_n_max()
    s_beam = s_beam.total_field(i_beam)

    a, b = i_beam.get_coefficients()
    p, q = s_beam.get_coefficients()
    n, m = i_beam.get_mode_indices()
    n_max = i_beam.get_n_max()
#    n = n.reshape((n.shape[0], 1)) if len(n.shape)==1 else m
#    m = m.reshape((m.shape[0], 1)) if len(m.shape)==1 else m
    b = 1j*b
    q = 1j*q
    addv = np.zeros((2*n_max+3,1))
    at = np.concatenate([a, matlib.repmat(addv, 1, a.shape[1])])
    bt = np.concatenate([b, matlib.repmat(addv, 1, b.shape[1])])
    pt = np.concatenate([p, matlib.repmat(addv, 1, p.shape[1])])
    qt = np.concatenate([q, matlib.repmat(addv, 1, q.shape[1])])
    ci = combined_index(n,m).astype(int)
    np1 = 2*n+1
    cinp1 = ci+np1
    cinp1mp1 = ci+np1+1
    cinp1mm1 = ci+np1-1
    cimp1 = ci+1
    kimp = (m>n-1).reshape((m.shape[0]))
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
    amp1 = at[ci, :]
    bmp1 = bt[ci, :]
    pmp1 = pt[ci, :]
    qmp1 = qt[ci, :]
    amp1[kimp, :] = 0
    bmp1[kimp, :] = 0
    pmp1[kimp, :] = 0
    qmp1[kimp, :] = 0 
    a = a[(ci-1), :]
    b=b[(ci-1), :]
    p=p[(ci-1), :]
    q=q[(ci-1), :]
    n = n.reshape((a.shape))# if len(n.shape)==1 else m
    m = m.reshape((a.shape))# if len(m.shape)==1 else m
    Az = m/n/(n+1)*(-(a)*b.conjugate()+(q.conjugate())*p).imag
    bz_part1 = 1/(n+1)*np.sqrt(n*(n-m+1)*(n+m+1)*(n+2)/(2*n+3)/(2*n+1))
    bz_part2 = (anp1*(a.conjugate())+bnp1*(b.conjugate())-(pnp1)*(p.conjugate())-qnp1*(q.conjugate())).imag
    Bz = bz_part1*bz_part2
    fz = 2*(Az.sum()+Bz.sum())
    Axy=1j/n/(n+1)*np.sqrt((n-m)*(n+m+1))*((pmp1.conjugate())*q - (amp1.conjugate())*b - (qmp1.conjugate())*p + (bmp1.conjugate())*a)
    bxy_part1 = 1j/(n+1)*np.sqrt(n*(n+2))/np.sqrt((2*n+1)*(2*n+3))
    bxy_part2 = np.sqrt((n+m+1)*(n+m+2))*(p*(pnp1mp1.conjugate()) + q*(qnp1mp1.conjugate()) -a*(anp1mp1.conjugate()) -b*(bnp1mp1.conjugate())) 
    bxy_part3 = np.sqrt((n-m+1)*(n-m+2))*(pnp1mm1*(p.conjugate()) + qnp1mm1*(q.conjugate()) - anp1mm1*(a.conjugate()) - bnp1mm1*(b.conjugate()))
    Bxy=bxy_part1*(bxy_part2+bxy_part3)
    fxy = (Axy+Bxy).sum()
    fx = fxy.real
    fy = fxy.imag
    return np.array([fx,fy,fz])