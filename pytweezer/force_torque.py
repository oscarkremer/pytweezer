import numpy as np
from .t_matrix import TMatrixMie
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
        n_beams = ibeam.n_beams
        if coherent:
            n_beams = 1
        f = np.zeros((3*T.T.shape, n_locations, nbeams))
        t = np.zeros((3*T.T.shape, n_locations, nbeams))
        s = np.zeros((3*T.T.shape, n_locations, nbeams))
        for i in range(1, n_locations+1):
            if position.size:
                if n_positions == 1:
                    aux_position = position
                else:
                    aux_position = position[:,i]
            if rotation.size:
                if n_rotations == 1:
                    aux_rotation = rotation
                else:
                    aux_rotation = rotation[:,i]
            print('computing scatter')
            sbeam, tbeam = ibeam.scatter(T, position=aux_position, rotation=aux_rotation)
            if coherent:
                sbeam = sbeam.mergeBeams()
                tbeam = tbeam.mergeBeams()
            fl, tl, sl = _force_torque_(tbeam, sbeam)
            f[:, i, :] = fl.reshape((3*T.size, 1, n_beams))
            t[:, i, :] = tl.reshape((3*T.size, 1, n_beams))
            z[:, i, :] = sl.reshape((3*T.size, 1, n_beams))
        f = squeeze(f)
        t = squeeze(t)
        s = squeeze(s)
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
'''  

def _force_torque_():
    % Check the number of beams in each input
    if ibeam.Nbeams ~= sbeam.Nbeams && ibeam.Nbeams ~= 1 && sbeam.Nbeams ~= 1
    error('Beam objects must contain same number of beams or 1 beam');
    end

    % Ensure beams are the same size
    if ibeam.Nmax > sbeam.Nmax
    sbeam.Nmax = ibeam.Nmax;
    elseif ibeam.Nmax < sbeam.Nmax
    ibeam.Nmax = sbeam.Nmax;
    end

    % Ensure the beam is incoming-outgoing
    sbeam = sbeam.totalField(ibeam);

    % Get the relevent beam coefficients
    [a, b] = ibeam.getCoefficients();
    [p, q] = sbeam.getCoefficients();
    [n, m] = ibeam.getModeIndices();

    nmax=ibeam.Nmax;

    b=1i*b;
    q=1i*q;

    addv=zeros(2*nmax+3,1);

    at=[a;repmat(addv, 1, size(a, 2))];
    bt=[b;repmat(addv, 1, size(b, 2))];
    pt=[p;repmat(addv, 1, size(p, 2))];
    qt=[q;repmat(addv, 1, size(q, 2))];

    ci=ott.utils.combined_index(n,m);

    %these preserve order and number of entries!
    np1=2*n+2;
    cinp1=ci+np1;
    cinp1mp1=ci+np1+1;
    cinp1mm1=ci+np1-1;
    cimp1=ci+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %this is for m+1... if m+1>n then we'll ignore!
    kimp=(m>n-1);

    anp1=at(cinp1, :);
    bnp1=bt(cinp1, :);
    pnp1=pt(cinp1, :);
    qnp1=qt(cinp1, :);

    anp1mp1=at(cinp1mp1, :);
    bnp1mp1=bt(cinp1mp1, :);
    pnp1mp1=pt(cinp1mp1, :);
    qnp1mp1=qt(cinp1mp1, :);

    anp1mm1=at(cinp1mm1, :);
    bnp1mm1=bt(cinp1mm1, :);
    pnp1mm1=pt(cinp1mm1, :);
    qnp1mm1=qt(cinp1mm1, :);

    amp1=at(cimp1, :);
    bmp1=bt(cimp1, :);
    pmp1=pt(cimp1, :);
    qmp1=qt(cimp1, :);

    amp1(kimp, :)=0;
    bmp1(kimp, :)=0;
    pmp1(kimp, :)=0;
    qmp1(kimp, :)=0;

    a=a(ci, :);
    b=b(ci, :);
    p=p(ci, :);
    q=q(ci, :);

    Az=m./n./(n+1).*imag(-(a).*conj(b)+conj(q).*(p)); %this has the correct sign... farsund. modes match.
    Bz=1./(n+1).*sqrt(n.*(n-m+1).*(n+m+1).*(n+2)./(2*n+3)./(2*n+1)) ... %.*n
        .*imag(anp1.*conj(a)+bnp1.*conj(b)-(pnp1).*conj(p) ...
        -(qnp1).*conj(q)); %this has the correct sign... farsund. modes match.

    fz=2*sum(Az+Bz);

    Axy=1i./n./(n+1).*sqrt((n-m).*(n+m+1)).*(conj(pmp1).*q - conj(amp1).*b - conj(qmp1).*p + conj(bmp1).*a); %this has the correct sign... farsund. modes match.
    Bxy=1i./(n+1).*sqrt(n.*(n+2))./sqrt((2*n+1).*(2*n+3)).* ... %sqrt(n.*)
        ( sqrt((n+m+1).*(n+m+2)) .* ( p.*conj(pnp1mp1) + q.* conj(qnp1mp1) -a.*conj(anp1mp1) -b.*conj(bnp1mp1)) + ... %this has the correct sign... farsund. modes match.
        sqrt((n-m+1).*(n-m+2)) .* (pnp1mm1.*conj(p) + qnp1mm1.*conj(q) - anp1mm1.*conj(a) - bnp1mm1.*conj(b)) ); %this has the correct sign... farsund. modes match.

    fxy=sum(Axy+Bxy);
    fx=real(fxy);
    fy=imag(fxy);

    if nargout > 1:
        tz=sum(m.*(a.*conj(a)+b.*conj(b)-p.*conj(p)-q.*conj(q))); %this has the correct sign... farsund. modes match.
        
        txy=sum(sqrt((n-m).*(n+m+1)).*(a.*conj(amp1)+b.*conj(bmp1)-p.*conj(pmp1)-q.*conj(qmp1))); %this has the correct sign... farsund. modes match.
        tx=real(txy);
        ty=imag(txy);
    
        if nargout > 2:
            Cz=m./n./(n+1).*(-(a).*conj(a)+conj(q).*(q)-(b).*conj(b)+conj(p).*(p));
            Dz=-2./(n+1).*sqrt(n.*(n-m+1).*(n+m+1).*(n+2)./(2*n+3)./(2*n+1)) ...
                .*real(anp1.*conj(b)-bnp1.*conj(a)-(pnp1).*conj(q) ...
                +(qnp1).*conj(p));
            
            sz = sum(Cz+Dz);
            
            Cxy=1i./n./(n+1).*sqrt((n-m).*(n+m+1)).*(conj(pmp1).*p - conj(amp1).*a + conj(qmp1).*q - conj(bmp1).*b);
            Dxy=1i./(n+1).*sqrt(n.*(n+2))./sqrt((2*n+1).*(2*n+3)).* ...
                ( (sqrt((n+m+1).*(n+m+2)) .* ( p.*conj(qnp1mp1) - q.* conj(pnp1mp1) -a.*conj(bnp1mp1) +b.*conj(anp1mp1))) + ...
                (sqrt((n-m+1).*(n-m+2)) .* (pnp1mm1.*conj(q) - qnp1mm1.*conj(p) - anp1mm1.*conj(b) + bnp1mm1.*conj(a))) );
            
            sxy=sum(Cxy+Dxy);
            sy=real(sxy);
            sx=imag(sxy);
    return fx,fy,fz,tx,ty,tz,sx,sy,sz
'''