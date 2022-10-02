'''
This script execute the kinematics and inverse kinematics computation.
Firstly, a set of joint angles must be inputed, used in the direct kinematic 
movement. Then, with the encountered position and orientation the inverse kinematics
is applied to verify if the same joint angles are found.
'''
from pymentor import Mentor
import numpy as np
import sys, os
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../')


def test_direct_kinematics():
    angles = [np.pi/6]*5
    robot = Mentor()
    pos, rot = robot.get_position(angles)
    print(pos, rot)
       
    
def test_inverse_kinematics1():
    pos = np.array([24.21323027, 13.97951501, -17.07885504, 1.0])
    rot = np.array([[0.59049287, 0.23642905, -0.77163428],
        [-0.23642905, -0.86349762, -0.44550326],
        [-0.77163428, 0.44550326, -0.4539905 ]])
    robot = Mentor()
    angles = robot.get_angles(pos,rot)
    pos, rot = robot.get_position(angles)
    print(pos,rot)

def test_inverse_kinematics2():
    robot = Mentor()
    pos = [10, 10, 10, 1]
    rot = robot.get_orientation(30, 30, 30)
    print(robot.get_angles(pos,rot)/np.pi)