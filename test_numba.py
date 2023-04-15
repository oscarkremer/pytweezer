from numba import jit, pycc

cc = pycc.CC('f')
@cc.export('', 'Tuple((i8, i8))(i8, i8)')
def f(x, y):
    # A somewhat trivial example
    return (x + y) / 3.14, x


f(1,1)
print(help(f))