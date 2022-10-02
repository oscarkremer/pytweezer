'''
This script defines the errors for the mentor class package.
The errors are related with impossibility to reach some positions and/or orientations.
'''

class InvalidOrientation(Exception):
    '''
    Exception when orientation entered is impossible to be reached
    by Mentor by mentor.
    '''
    pass


class InvalidPosition(Exception):
    '''
    Exception to be raised when entered position is impossible to be reached
    by Mentor.
    '''
    pass


class InvalidPair(Exception):
    '''
    Exception to be raised when position and orientation entered are impossible to be reached
    by Mentor.
    '''
    pass