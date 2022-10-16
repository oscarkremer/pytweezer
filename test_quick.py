from pytweezer.utils import translate_z
import numpy as np
A, B, C = translate_z(1, -np.pi/2, method='videen', function_type='sbesselh2')


print(A, B, C)

print('-----')
for i in range(C.shape[2]):
    print(C[:,:,i])