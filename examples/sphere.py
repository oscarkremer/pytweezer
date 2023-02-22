from pytweezer.t_matrix import TMatrixMie
from pytweezer.beams import Gaussian
from pytweezer import force_torque
n_medium = 1.33
n_particle = 1.59

wavelength0 = 1064e-9

wavelength_medium = wavelength0/n_medium

radius = 1.0*wavelength_medium

beam_type = 'gaussian';

NA = 1.02;

T = ott.Tmatrix.simple('sphere', radius, 'wavelength0', wavelength0, ...
    'index_medium', n_medium, 'index_particle', n_particle);


beam = ott.BscPmGauss('NA', NA, 'polarisation', [ 1 1i ], ...
        'index_medium', n_medium, 'wavelength0', wavelength0);

beam.power = 1.0;

disp(['Beam calculation took ' num2str(toc) ' seconds']);

z = [0;0;1]*linspace(-8,8,80)*wavelength_medium;
fz = ott.forcetorque(beam, T, 'position', z);

zeq = ott.find_equilibrium(z(3, :), fz(3, :));
if isempty(zeq)
  warning('No axial equilibrium in range!')
  zeq=0;
zeq = zeq(1);

r = [1;0;0]*linspace(-4,4,80)*wavelength_medium + [0;0;zeq];
fr = ott.forcetorque(beam, T, 'position', r);

disp(['Force calculation took ' num2str(toc) ' seconds']);

%% Generate the plots

figure(1); plot(z(3, :)/wavelength_medium,fz(3, :));
xlabel('{\it z} [\lambda_m]');
ylabel('{\it Q_z} [n_m P / c]');
aa = axis;
hold on;
line(aa(1:2),[ 0 0 ],'linestyle',':');
line([0 0],aa(3:4),'linestyle',':');
hold off;

figure(2); plot(r(1, :)/wavelength_medium,fr(1, :));
xlabel('{\it r} [\lambda_m]');
ylabel('{\it Q_r} [n_m P / c]');
aa = axis;
hold on;
line(aa(1:2),[ 0 0 ],'linestyle',':');
line([0 0],aa(3:4),'linestyle',':');
hold off;
