import sys; sys.path.insert(0,'.')
from rcore import *
B=poly_of([(0,1),(0,2)])*12345678901234567890123-poly_of([(1,2),(3,4)])*98765432109876543210
c,f=factorZ(B); print(c,f)
import sympy
print(full_gain([(0,1)]*2,[(0,2),(1,2)],-1,4))
