# Affine dimension of the convex hull of the 20 code vertices x_R, iota x_R (|R|=3) inside Pi_nt (dim 11),
# and whether 1/2 is a relative interior point; also add the even-set code vertex F_[5] and its complement.
import itertools, numpy as np, sympy as sp
from scipy.optimize import linprog
Bs=[frozenset(c) for r in range(1,5) for c in itertools.combinations(range(5),r)]
def code(R,parity):
    return np.array([1.0 if len(B&set(R))%2==parity else 0.0 for B in Bs])
pts=[]
for R in itertools.combinations(range(5),3):
    pts.append(code(R,0)); pts.append(code(R,1))
P=np.array(pts)
D=sp.Matrix(P[1:]-P[0])
print("affine dim of 20 code vertices (|R|=3):",D.rank())
P2=np.vstack([P,code(range(5),0),code(range(5),1)])
print("affine dim with F_[5] and complement:",sp.Matrix(P2[1:]-P2[0]).rank())
# is 1/2 in the relative interior of conv(P)? check that every direction in the affine hull can be moved: LP max/min
half=np.full(30,0.5)
# write 1/2 as convex combination with all weights >= eps (strictly positive) -> relative interior
n=len(P)
res=linprog(-np.ones(1),A_ub=None,b_ub=None,A_eq=None,b_eq=None,bounds=[(None,None)]) if False else None
# maximize t s.t. sum l_i P_i = 1/2, sum l_i = 1, l_i >= t
c=np.zeros(n+1); c[-1]=-1
A_eq=np.vstack([np.hstack([P.T,np.zeros((30,1))]),np.hstack([np.ones((1,n)),np.zeros((1,1))])])
b_eq=np.concatenate([half,[1]])
A_ub=np.hstack([-np.eye(n),np.ones((n,1))]); b_ub=np.zeros(n)
r=linprog(c,A_ub=A_ub,b_ub=b_ub,A_eq=A_eq,b_eq=b_eq,bounds=[(None,None)]*(n+1),method='highs')
print("max min-weight of 1/2 as convex combination of the 20 code vertices:",r.x[-1] if r.status==0 else r.status)
