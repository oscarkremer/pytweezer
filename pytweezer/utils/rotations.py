import numpy as np


def rotx(angle_deg):
  '''
  Create a 3x3 rotation matrix for rotation about x axis
  
  R = rotx(angle_deg) calculate the rotation matrix for rotations from
  the +z towards +y axis.
  '''
  theta = angle_deg * np.pi/180
  return np.array([[1, 0, 0],
          [0, np.cos(theta), -np.sin(theta)],
          [0, np.sin(theta), np.cos(theta)]])

def roty(angle_deg):
  '''
  Create a 3x3 rotation matrix for rotation about x axis
  
  R = roty(angle_deg) calculate the rotation matrix for rotations from
  the +z towards +x axis.
  '''
  theta = angle_deg * np.pi/180
  return np.array([[np.cos(theta), 0, np.sin(theta)],
          [0, 1, 0],
          [-np.sin(theta), 0, np.cos(theta)]])

def rotz(angle_deg):
  '''
  Create a 3x3 rotation matrix for rotation about x axis
  
  R = rotx(angle_deg) calculate the rotation matrix for rotations from
  the +x towards +y axis.
  '''
  theta = angle_deg * np.pi/180
  return np.array([[np.cos(theta), -np.sin(theta), 0], 
                  [np.sin(theta), np.cos(theta), 0],
                  [0, 0, 1]])