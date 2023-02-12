from pytweezer.shapes import Sphere
import numpy as np

radius = 1.0
offset = np.array([0, 0, 2])
shape = Sphere(radius, position=offset)
z = np.linspace(-10, 10, 100)
mask = shape.inside_xyz(0, 0, z.T, origin='world') 
