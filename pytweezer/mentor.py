'''
This script defines the Mentor class package, where there are
specified as Mentor methods the computations of direct and
inverse kinematics
'''
import itertools
import numpy as np
from pymentor.error import InvalidOrientation, InvalidPosition, InvalidPair

TOL = 0.001

class Mentor:
    '''
    Mentor Class. This class represents the mentor didactic robot, obeying 
    the laws of direct and inverse kinematics, and any physical constraint
    of movement, both angulars position or velocity. 

    Attributes
    ----------
    alpha : list
        Angle list of Denavit Hartenberg parameter.
    a: list
        List of perpendicular distance between joints.
    d: list
        List of perpendicular distance between links.
    
    Methods
    -------
    test_inverse_kinematics(self, pos, rot, tag_theta1=True, tag_theta2=True, tag_theta3=True):  
        Method to test with a certain set of condition represents a possible situation for the robot.
    get_angles(self, pos, rot):
        Method to test with all possible scenarious of angles.
    _inverse_kinematics(self, pos, orientation, tag_theta1=True, tag_theta2=True, tag_theta3=True):
        Angles computation in inverse kinematics.
    verify(self, pos, rot, returned_pos, returned_rot):
        Verification if computed position is equal to desired position.
    get_position(self, theta, z_axis=0):
        Method to apply direct kinematics 
    get_orientation(self, alpha, beta, gamma):
        Method to get end effector orientation consideration Alpha, Beta and
        Gamma - XYZ angles.
    separate(self, matrix, z_axis):
        Method to separate a 4x4 matrix in rotational matrix and position array.
    denavit(self, theta, i):
        Method to apply the Denavit-Hartenberg to find the transformation matrix.
    fix_quadrante(self, sin, cos):
        Method to adjust the angle quadrant from its sine and cosine.
    fix_theta(self, param, tag, angle):
        Method to adjust the angle quadrant from any possible set o parameters.
    '''
    def __init__(self):
        '''
        The constructor method doesn't use any parameters. 
        This method is only responsible to define the Denavit-Hartenberg 
        parameters as Mentor attributes.
        '''
        self.alpha = [0, -np.pi/2, 0, 0, -np.pi/2]
        self.a = [0, 0, 17.2739, 15.5, 0]
        self.d = [0, 0, 0, 0, 0]

    def test_inverse_kinematics(self, pos, rot, tag_theta1=True, tag_theta2=True, tag_theta3=True):  
        '''
        Method to test with a certain set of condition represents a possible situation for the robot.

        Parameters
        ----------
        pos: list
            List with position in the cartesian space of the end effector grip.
        rot: list
            Rotation matrix of the end effector grip.
        tag_theta1: bool(optional):
            Tag to identify error in theta-1
        tag_theta2: bool(optional):
            Tag to identify error in theta-2
        tag_theta3: bool(optional):
            Tag to identify error in theta-3
 
        Returns
        -------
        tag:
            Verification if there was any error in inverse kinematics.
        theta:
            List containing the joint angles encountered.
        '''

        theta = self._inverse_kinematics(pos, rot, tag_theta1=tag_theta1, tag_theta2=tag_theta2, tag_theta3=tag_theta3)
        returned_pos, returned_rot = self.get_position(theta)        
        tag_orientation, tag_position, tag_pair = self.verify(pos, rot, returned_pos, returned_rot)
        return tag_orientation, tag_position, tag_pair, theta

    def get_angles(self, pos, rot):
        '''
        Method to test with all possible scenarious of angles.

        Parameters
        ----------
        pos: list
            List with position in the cartesian space of the end effector grip.
        rot: list
            Rotation matrix of the end effector grip.

        Returns
        -------
        tag:
            Tag that identifies if position is or isn't possible.
        theta:
            List of joint angles returned in case position/orientation is possible.
        '''
        for element in itertools.product([True, False],[False, True],[False, True]):
            tag_orientation, tag_position, tag_pair, theta = self.test_inverse_kinematics(pos, rot, element[0], element[1], element[2])
            if tag_orientation*tag_position*tag_pair:
                return theta
        if tag_pair:
            raise InvalidPair('Pair Position-Orientation unreachable, please try different values!')
        if tag_orientation:
            raise InvalidOrientation('Orientation unreachable for such position, please try different values for alpha-beta-gamma angles!')
        if tag_position:
            raise InvalidPosition('Position unreachable, please try different values for x-y-z coordinates!')

    def _inverse_kinematics(self, pos, orientation, tag_theta1=True, tag_theta2=True, tag_theta3=True):
        '''
        Angles computation in inverse kinematics.

        Parameters
        ----------
        pos: list
            List with position in the cartesian space of the end effector grip.
        rot: list
            Rotation matrix of the end effector grip.
        tag_theta1: bool(optional):
            Tag to identify error in theta-1.
        tag_theta2: bool(optional):
            Tag to identify error in theta-2.
        tag_theta3: bool(optional):
            Tag to identify error in theta-3.
         
        Returns
        -------
        theta:
            List of joint angles returned in case position/orientation is possible.
        '''
        theta = []
        theta1 = np.nan_to_num(self.fix_theta(pos, tag_theta1, 'theta1'))
        theta.append(theta1)
        theta3 = self.fix_theta((pos[0]**2+pos[1]**2+pos[2]**2-self.a[2]**2-self.a[3]**2)/(2*self.a[2]*self.a[3]), tag_theta3, 'theta3')
        theta3 = np.nan_to_num(theta3)
        theta2  = np.nan_to_num(-theta3+self.fix_theta(((-np.cos(theta3)*self.a[2]-self.a[3])*pos[2]+np.sin(theta3)*self.a[2]*(np.sin(theta1)*pos[1]+pos[0]*np.cos(theta1)))/(np.power(np.cos(theta1)*pos[0]+np.sin(theta1)*pos[1], 2)+pos[2]*pos[2]), tag_theta2, 'theta2'))
        theta.append(theta2)
        theta.append(theta3)
        sin_theta4 = np.sin(theta2+theta3)*orientation[2][2] - np.cos(theta1)*np.cos(theta2+theta3)*orientation[0][2] - np.sin(theta1)*np.cos(theta2+theta3)*orientation[1][2]
        cos_theta4 = -np.cos(theta2+theta3)*orientation[2][2] - np.cos(theta1)*np.sin(theta2+theta3)*orientation[0][2] - np.sin(theta1)*np.sin(theta2+theta3)*orientation[1][2]
        theta.append(self.fix_quadrante(sin_theta4, cos_theta4))
        sin_theta5 = np.sin(theta1)*orientation[0][0] - np.cos(theta1)*orientation[1][0]
        cos_theta5 = np.sin(theta1)*orientation[0][1] - np.cos(theta1)*orientation[1][1]
        theta.append(self.fix_quadrante(sin_theta5, cos_theta5))
        return theta

    def verify(self, pos, rot, returned_pos, returned_rot):
        '''
        This method executes verification if computed position is equal to the desired position.

        Parameters
        ----------
        pos: list
            List with position in the cartesian space of the end effector grip.
        rot: list
            Rotation matrix of the end effector grip.
        returned_pos: list
            List with returned position if ones applies the angles encontered 
            by the inverse kinematics in the direct kinematics.
        returned_rot: np.array
            3x3 matrix representing rotation encountered if ones applies the angles encontered 
            by the inverse kinematics in the direct kinematics.
        
        Returns
        -------
        Tag identifying if rotation and position returned matches with the ones pre-specified.
        '''
        orientation, position, pair = False, False, False
        diff = pos[:3]-returned_pos[:3]
        rot_norm = np.linalg.norm(rot-returned_rot)
        if rot_norm > TOL and (abs(diff[0])>TOL or abs(diff[1])>TOL or abs(diff[2])>TOL):
            return False, False, False
        else:
            if rot_norm > TOL:
                return False, True , True            
            else:
                if (abs(diff[0])>TOL or abs(diff[1])>TOL or abs(diff[2])>TOL):
                    return True, False, True
                else:
                    return True, True, True



    def get_position(self, theta, z_axis=0):
        '''
        Method to apply direct kinematics in a set of inputed joint angles.

        Parameters
        ----------
        theta: Numpy Array
            Array with joint angles of the robot.
        z_axis: float (optional)
            Orthogonal distance to create an extra axis located with the same orientation 
            of the frame from the theta 5 but with a displacement in the z-axis.
            
        Returns
        -------
        Position array and rotation matrix which define the end effector grip in 
        the cartesian space.
        '''
        matrix = np.matmul(self.denavit(theta, 3),self.denavit(theta, 4))
        for i in range(3):
            matrix = np.matmul(self.denavit(theta, 2-i), matrix)
        return self.separate(matrix, z_axis)
    
    def get_orientation(self, alpha, beta, gamma):
        '''
        Method to get end effector orientation consideration Alpha, Beta and
        Gamma - XYZ angles.

        Parameters
        ----------
        alpha: float
            Orientation angle of X axis.
        beta: float
            Orientation angle of Y axis.
        gamma: float
            Orientation angle of Z axis.

        Returns
        -------
        Rotation 3x3 matrix created from the alpha, beta and gamma angles specified in 
        the input functions.
        '''
        return [[np.cos(alpha)*np.cos(beta), np.cos(alpha)*np.sin(beta)*np.sin(gamma) - np.sin(alpha)*np.cos(gamma), np.cos(alpha)*np.sin(beta)*np.cos(gamma) + np.sin(alpha)*np.sin(gamma) ],
                        [np.sin(alpha)*np.cos(beta), np.sin(alpha)*np.sin(beta)*np.sin(gamma) + np.cos(alpha)*np.cos(gamma), np.sin(alpha)*np.sin(beta)*np.cos(gamma) - np.cos(alpha)*np.sin(gamma) ],
                        [-np.sin(beta), -np.cos(beta)*np.sin(gamma), np.cos(beta)*np.cos(gamma)]]

    def separate(self, matrix, z_axis):
        '''
        Method to separate a 4x4 matrix in rotational matrix and position array.

        Parameters
        ----------
        Matrix: Numpy Array
            Joint angle 4x4 matrix.
        z_axis: float (optional)
            Distance in Z Axis to define an extra frame.

        Returns
        -------
        pos: float
            3x1 Array that defines the position of end effector grip.
        np.array(rot):
            3x3 Matrix that defines the orientation of end effector grip.
        '''
        pos =  np.matmul(matrix, [0,0,z_axis,1])
        rot = [matrix[i][:-1] for i in range(3)]
        return pos, np.array(rot)

    def denavit(self, theta, i):
        '''
        Method to apply the Denavit-Hartenberg to find the transformation matrix.

        Parameters
        ----------
        theta: list
            List of joint angles.    
        i: int
            Index to determine which transformation matrix is being calculated.
            The index i represents transformation from i to i+1 frame.

        Returns
        -------
        4x4 Block Triangular Matrix which contain both Rotation matrix (3x3) and 
            the position vector (3x1).
        '''
        return  [[np.cos(theta[i]), -np.sin(theta[i]),  0,   self.a[i]],
            [np.sin(theta[i])*np.cos(self.alpha[i]), np.cos(theta[i])*np.cos(self.alpha[i]), -np.sin(self.alpha[i]), -np.sin(self.alpha[i])*self.d[i]],
            [np.sin(theta[i])*np.sin(self.alpha[i]), np.cos(theta[i])*np.sin(self.alpha[i]), np.cos(self.alpha[i]), np.cos(self.alpha[i])*self.d[i]],
            [0, 0, 0, 1]]

    def fix_quadrante(self, sin, cos):
        '''
        Method to adjust the angle quadrant from its sine and cosine.

        Parameters
        ----------
        sin: float
            Angle sine.
        cos: float
            Angle cossine.
        
        Returns
        -------
        Angle adjusted to the right quadrant.
        '''
        if sin >= 0 and cos >= 0:
            return np.arcsin(abs(sin))
        if sin >= 0  and cos <= 0:
            return np.pi - np.arcsin(abs(sin))
        if sin <= 0 and cos >= 0:
            return 2*np.pi - np.arcsin(abs(sin))
        if sin <= 0 and cos <= 0:
            return np.arccos(abs(cos)) + np.pi

    def fix_theta(self, param, tag, angle):
        '''
        Method to adjust the angle quadrant in a more generic way.

        Parameters
        ----------
        param: float
            Used parameter to identify if the angle is correct. Param can be
            a list of (x, y), sine or cossine.
        tag: boolean
            Tag that identifies if any error have already been detected.
        angle: str
            Specific angle analysed.
            
        Returns
        -------
        Angle adjusted to the right quadrant.
        '''
        if angle=='theta1':
            if abs(param[0]) < 0.001:  
                param = np.pi/2
            else:
                param = np.arctan(param[1]/ param[0])
        if tag and param>=0:
            if angle == 'theta1':
                return param
            if angle == 'theta2':
                return np.arcsin(abs(param))       
            if angle == 'theta3':
                return np.arccos(abs(param))       
        if not tag and param>=0:
            if angle == 'theta1':
                return np.pi + param
            if angle == 'theta2':
                return np.pi-np.arcsin(abs(param))
            if angle == 'theta3':
                return 2*np.pi-np.arccos(abs(param))   
        if tag and param<=0:
            if angle == 'theta1':
                return 2*np.pi - abs(param)
            if angle == 'theta2':
                return 2*np.pi - np.arcsin(abs(param))
            if angle == 'theta3':   
                return np.pi-np.arccos(abs(param))       
        if not tag and param<=0:
            if angle == 'theta1':
                return np.pi - abs(param)
            if angle == 'theta2':
                return np.pi + np.arcsin(abs(param))
            if angle == 'theta3':
                return np.pi + np.arccos(abs(param))