# Walsh block profile of VI.5 (one column per block; each block's primes share the column)
from sympy import Matrix, lcm
from sympy.matrices.normalforms import hermite_normal_form
import itertools
for k in range(2,6):
    m=2**k
    vs=[v for v in itertools.product([0,1],repeat=k) if any(v)]
    rows=[[sum(a*b for a,b in zip(l,v))%2 for v in vs] for l in itertools.product([0,1],repeat=k)]
    B=Matrix([[2*(r[j]-rows[0][j]) for j in range(m-1)] for r in rows[1:]])
    Bas=hermite_normal_form(B.T).T
    Bas=Matrix([list(Bas.row(i)) for i in range(Bas.rows) if any(Bas.row(i))])
    inv=Bas.inv()
    ords=set()
    for j in range(m-1):
        e=Matrix([[1 if i==j else 0 for i in range(m-1)]]); x=e*inv
        ords.add(int(lcm([xx.q for xx in x])))
    print("M=",m,"rank",B.rank(),"index",abs(Bas.det()),"orders of block indicators",ords)
