"""Two nonzero Jacobian minors certify all 16 two-moment points."""
from itertools import combinations
from sympy import Matrix, Symbol, rem
w=Symbol('w')
edges=list(combinations(range(5),2))
A=Matrix([[int(i in e) for e in edges] for i in range(5)])
cycle=[1 if (j-i)%5 in (1,4) else -1 for i,j in edges]
tetra=[0 if 4 in (i,j) else 1 if (i,j) in ((0,1),(2,3)) else w if (i,j) in ((0,2),(1,3)) else -1-w for i,j in edges]
def minor(values,cols):
    B=Matrix([[(int(i in e)-int(4 in e))*y for e,y in zip(edges,values)] for i in range(4)])
    return A.col_join(B)[:,cols].det().expand()
cycleminor=minor(cycle,list(range(9)))
tetminor=rem(minor(tetra,[0,1,2,3,4,5,6,8,9]),w*w+w+1,w)
print('cycle minor',cycleminor)
print('tetrahedral minor',tetminor)
