import numpy as np
def xyz2rtp(*args):
   '''
   % XYZ2RTP coordinate transformation from cartesian to spherical
   %   r      radial distance [0, Inf)
   %   theta  polar angle, measured from +z axis [0, pi]
   %   phi    azimuthal angle, measured from +x towards +y axes [0, 2*pi)
   %
   % [r,theta,phi] = XYZ2RTP(x,y,z) takes vectors or scalars outputs
   % the spherical coordinates as vectors/scalars of the same size.
   %
   % [r,theta,phi] = XYZ2RTP(x) same as above but with the coordinate
   % packed into the vector/matrix x = [ x y z ].
   %
   % r = XYZ2RTP(...) same as above with the result packed into
   % the vector/matrix r = [ r theta phi ].

   % This file is part of the optical tweezers toolbox.
   '''

   def __mod__(a, m):
      return a - m*np.floor(a/m)

   if len(args) == 1:
      array_input = args[0]
      y = array_input[1]
      z = array_input[2]
      x = array_input[0]
   else:
      x = args[0]
      y = args[1]
      z = args[2]

   xy = np.sqrt(np.power(x, 2) + np.power(y, 2))
   theta = __mod__(np.arctan2(xy, z)+2*np.pi, 2*np.pi)
   phi = __mod__(np.arctan2(y, x)+2*np.pi, 2*np.pi)
   r = np.sqrt(x*x + y*y + z*z)
   return r.T.reshape(-1), theta.T.reshape(-1), phi.T.reshape(-1)
