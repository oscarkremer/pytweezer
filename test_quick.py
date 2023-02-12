from pytweezer.shapes import Sphere
import numpy as np

decimal = 3
radius = 1 
numr = 100
sphere_1 = Sphere(radius)
sphere_2 = Sphere(radius, position=np.array([1,0,0]))
spacing = 0.1
voxels1 = sphere_1.voxels(spacing, origin='world') + offset
voxels2 = sphere_2.voxels(spacing, origin='world')