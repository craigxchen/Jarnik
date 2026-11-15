import numpy as np, itertools, flint
P=np.array([[1,1,1,1,1,-1,-1,-1],[1,1,1,-1,-1,1,1,-1],[1,1,-1,1,-1,1,-1,1],[1,-1,1,1,-1,-1,1,1],[1,-1,-1,-1,1,1,1,1],[-1,1,1,-1,1,1,-1,1]])
D=8;L=3
d=[(P[i]-P[0])//2 for i in range(1,6)]
M=flint.fmpz_mat([[int(x) for x in v] for v in d])
Y,nul=M.nullspace()
Yt=flint.fmpz_mat([[Y[r,c] for r in range(D)] for c in range(nul)])  # nul x D
Z,k=Yt.nullspace()  # D x k basis of V cap Q^D (integer, maybe not saturated)
B=flint.fmpz_mat([[Z[r,c] for r in range(D)] for c in range(k)])
# saturate: Hermite/LLL on the lattice V cap Z^D: use that V cap Z^D = kernel lattice of Yt; flint nullspace returns a Q-basis; saturate via HNF trick
import sympy
Ys=sympy.Matrix(nul,D,lambda i,j:int(Yt[i,j]))
# integer kernel basis via Smith-like: use sympy's nullspace + LLL on scaled lattice
from sympy.matrices.normalforms import hermite_normal_form
# Lattice {a in Z^D: Ys a = 0}: compute via LLL on [I | W*Ys^T]
W=10**6
Lat=flint.fmpz_mat([[ (1 if i==j else 0) for j in range(D)]+[W*int(Ys[r,i]) for r in range(nul)] for i in range(D)])
R=Lat.lll()
basis=[]
for i in range(D):
    row=[int(R[i,j]) for j in range(D+nul)]
    if all(x==0 for x in row[D:]): basis.append(row[:D])
basis=np.array(basis); print("saturated basis rank",len(basis)); print(basis)
mins=None;kill=None
for co in itertools.product(range(-4,5),repeat=len(basis)):
    a=np.array(co)@basis
    if not a.any(): continue
    ap=a[a>0].sum(); am=-a[a<0].sum()
    if ap+am<=L or 2*min(ap,am)<L: kill=a; break
print("kill",kill)
# index of span(d) in saturation
print("index", abs(sympy.Matrix(basis).T.solve_least_squares(sympy.Matrix(np.array(d).T)).det()) if len(basis)==5 else None)
