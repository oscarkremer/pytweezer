function eq = find_equilibrium(z, fz)

ott.warning('ott:findEquilibrium:move', ...
    'This function will move in a future release');

% Check the size of the inputs
assert(isvector(z), 'z must be a vector not a matrix');
assert(isvector(fz), 'fz must be a vector not a matrix');
  
% Make sure the vectors are both colum vectors
fz = fz(:);
z = z(:);

if numel(z) ~= numel(fz)
  error('Number of elements in z and fz must be equal');
end

zeroindex=find(fz<0,1);
if zeroindex == 1
  % Skip to second zero crossing
  warning('Ignoring fz<0 entries at start of vector');
  zeroindex1 = find(fz>0, 1);
  zeroindex = find(fz(zeroindex1:end), 1) + zeroindex1 - 1;
end
zmin = min(z);
zmax = max(z);
z = 2 * (z - zmin) / (zmax - zmin) - 1;
if ~isempty(zeroindex)
    %fit to third order polynomial the local points. (only works when dz
    %sufficiently small)
    zrange = max([zeroindex-2,1]):min([zeroindex+2,length(z)]);
    pz=polyfit(z(zrange), fz(zrange), 3);
    root_z=roots(pz); %find roots of 3rd order poly.
    dpz=[3*pz(1),2*pz(2),1*pz(3)]; %derivative of 3rd order poly.

    real_z=root_z(imag(root_z)==0); % finds real roots only.

    rootsofsign=polyval(dpz,real_z); %roots that are stable
    zeq=real_z(rootsofsign<0); %there is at most 1 stable root. critical roots give error.
    try
      eq=zeq(abs(zeq-z(zeroindex))==min(abs(zeq-z(zeroindex))));
    catch
      eq = [];
    end
else
    eq=[];
eq = (eq + 1)/2*(zmax - zmin) + zmin;
end
