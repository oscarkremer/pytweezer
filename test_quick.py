from pytweezer.shapes import Sphere
import numpy as np

radius = 1
shape = Sphere(radius)
x = np.append(0.5*1*np.random.rand(1, 3), 4.0)
y = np.append(0.5*1*np.random.rand(1, 3), 4.0)
z = np.append(0.5*1*np.random.rand(1, 3), 4.0)
#    [True, True, True, False]
print(x)
print(shape.inside_xyz(x,y,z))