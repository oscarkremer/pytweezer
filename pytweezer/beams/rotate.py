
from pytweezer.utils.rotation import *
from pytweezer.utils import combined_index, wigner_rotation_matrix

def rotate_yz(beam, angley, anglez, n_max)
    beam, D = rotate(np.matmul(rotz(anglez*180/np.pi)*roty(angley*180/np.pi)))
    return beam, D


def rotate(beam, wigner=np.array([]), R=np.array([]), n_max=0)
    n_max = beam.get_n_max() if not n_max else n_max
    if R.size and not wigner.size:
        if sum(sum((eye(3) - R).^2)) < 1e-6:
          D = eye(combined_index(p.Results.Nmax, p.Results.Nmax));
          return;
        D = ott.utils.wigner_rotation_matrix(...
            max(beam.Nmax, p.Results.Nmax), R);
        beam = rotate('wigner', D);
    elif ~isempty(p.Results.wigner) && isempty(p.Results.R):
        
        if iscell(p.Results.wigner):
            ibeam = beam;
            beam = ott.Bsc();

            for ii = 1:numel(p.Results.wigner):
                sz = size(ibeam.a, 1);
                D2 = p.Results.wigner{ii}(1:sz, 1:sz);
                beam = beam.append(D2 * ibeam);        
        else:
            sz = size(beam.a, 1);
            D2 = p.Results.wigner(1:sz, 1:sz);
            beam = D2 * beam;
      else:
        raise ValueError('One of wigner or R must be specified')
