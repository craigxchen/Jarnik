"""The Jacobian of the 4+1 singular degree-one branch at 7."""
from itertools import combinations
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix
from sympy import GF

p=7
edges=list(combinations(range(5),2))
a,b=(0,1),(2,6)
u,v=(4,4),(5,3)
def bar(z):
    return z[0],(-z[1])%p
roots=(bar(a),u,v,a,v,u,a,bar(b),b,b)
M=[]
for k in range(1,6):
    for i in range(4):
        row=[]
        for e,(x,y) in zip(edges,roots):
            A,B=1,0
            for _ in range(k-1):A,B=(A*x-B*y)%p,(A*y+B*x)%p
            s=int(i in e)-int(4 in e)
            row += [s*k*B,s*k*A]
        M.append(row)
def matrix(rows):
    return DomainMatrix.from_list(rows,ZZ).convert_to(GF(p))
mat=matrix(M)
print('rank',mat.rank())
N=mat.nullspace().to_Matrix()
perm={0:1,1:0,2:3,3:2,4:4}
pidx=[edges.index(tuple(sorted(perm[j] for j in e))) for e in edges]
print('residue',list(zip(edges,roots)))
print('simultaneous twin swap', all(roots[j]==roots[k] for j,k in enumerate(pidx)))
for j,e in enumerate(edges):
    h=pidx[j]
    print('fixed differential',e, N[:,2*j]==N[:,2*h] and N[:,2*j+1]==N[:,2*h+1])
