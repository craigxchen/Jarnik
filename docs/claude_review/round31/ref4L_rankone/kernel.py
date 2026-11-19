import sympy as sp
E=[(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]
M=sp.Matrix([[1 if i in e else 0 for e in E] for i in range(5)])
from sympy.matrices.normalforms import smith_normal_form
print(smith_normal_form(M,domain=sp.ZZ))
ns=M.nullspace()
B=sp.Matrix.hstack(*ns); B=B*sp.ilcm(*[x.q for x in B])
print(B.T)
