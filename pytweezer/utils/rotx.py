import numpy as np


def rotx(angle_deg):
  '''
  % Create a 3x3 rotation matrix for rotation about x axis
  %
  % R = rotx(angle_deg) calculate the rotation matrix for rotations from
  % the +z towards +y axis.
  %
  % R = rotx([a1, a2, a3, ...]) returns a 3xN matrix of rotation matrices
  % for each angle in the input.
  %
  % Optional named arguments:
  %   usecell    bool     True to output as cell array instead of 3xN matrix.
  %       Default: false.  The cell array has the same shape as angle_deg.
  %
  % Replacement/extension to Matlab rotx function provided in the
  % Phased Array System Toolbox.

  % This file is part of the optical tweezers toolbox.
  % See LICENSE.md for information about using/distributing this file.
  '''
  theta = angle_deg * np.pi/180;
  return np.array([[1, 0, 0],
          [0, np.cos(theta), -np.sin(theta)],
          [0, np.sin(theta), np.cos(theta)]])