import numpy as np
from stl import mesh
from pytweezer.shapes import STLLoader

cube = STLLoader('tests/shapes/cube.stl')
print(cube._verts)