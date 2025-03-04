from pytweezer.beams import Gaussian
from pytweezer.t_matrix import TMatrixMie
from pytweezer import force_torque, find_equilibrium
import time
import numpy as np
import matplotlib.pyplot as plt
n_medium = 1.33
n_particle = 1.59
wavelength0 = 1064e-9
wavelength_medium = wavelength0/n_medium
radius = 1.0*wavelength_medium
beam_type = 'gaussian'
NA = 1.02

T = TMatrixMie(radius=radius, lambda_0=wavelength0, index_m=n_medium, index_p=n_particle)

beam = Gaussian(power=1.0, na=NA, polarization=np.array([1, 1j]), index_m=n_medium, lambda_0=wavelength0)
start = time.time()#process_time()
#print(np.abs(fr).sum(), np.abs(fr).mean())
z = np.array([[0],[0],[1]])*np.linspace(-8,8,200)*wavelength_medium
f = force_torque(beam, T, position=z)

zeq = find_equilibrium(z[2,:], f[2, 0, :])[0].real
r = np.array([[1],[0],[0]])*np.linspace(-4,4,200)*wavelength_medium + np.array([[0],[0],[zeq]])

fr = force_torque(beam, T, position=r)
plt.plot(r[0,:]/wavelength_medium, fr[0,0,:])
plt.xlim([-8, 8])
plt.show()
