# Max size of a binary code C = V cap {0,1}^D, V a real LINEAR subspace of dim <= r, with min distance >= D/2 (asymptotic Plotkin regime).
# Codewords <-> S subset {0,1}^r (projection to an information set), coordinates <-> linear functionals l with l(S) in {0,1}.
import itertools, numpy as np
from fractions import Fraction
from scipy.optimize import linprog
import sys
r=int(sys.argv[1])
pts=list(itertools.product([0,1],repeat=r))
def functionals(S):
    S=np.array(S,dtype=float)
    # candidate l: solve l on basis subsets of S; simpler: l in Z^r with small entries? use all l with values determined by 0/1 on unit vectors? not general.
    # general: pick r linearly independent points B of S, assign 0/1 values, solve.
    out=set()
    idx=list(range(len(S)))
    for B in itertools.combinations(idx,r):
        MB=S[list(B)]
        if abs(np.linalg.det(MB))<1e-9: continue
        for vals in itertools.product([0,1],repeat=r):
            l=np.linalg.solve(MB,np.array(vals,float))
            v=S@l
            if np.all(np.abs(v-np.round(v))<1e-9) and np.all((np.round(v)==0)|(np.round(v)==1)):
                out.add(tuple(np.round(v).astype(int)))
        break  # one basis suffices: l determined by values on a basis
    return [o for o in out if any(o)]
def best_t(S):
    F=functionals(S)
    M=len(S); pairs=list(itertools.combinations(range(M),2))
    if not F: return 0
    # maximize t: sum_f w_f [f_i!=f_j] >= t, sum w=1
    A=[];b=[]
    for (i,j) in pairs:
        A.append([-int(f[i]!=f[j]) for f in F]+[1]); b.append(0)
    c=[0]*len(F)+[-1]
    res=linprog(c,A_ub=A,b_ub=b,A_eq=[[1]*len(F)+[0]],b_eq=[1],bounds=[(0,None)]*len(F)+[(None,None)])
    return -res.fun
best={}
rest=pts[1:]
for M in range(16,1,-1):
  found=False
  for comb in itertools.combinations(rest,M-1):
    S=[pts[0]]+list(comb)
    if np.linalg.matrix_rank(np.array(S))<r: continue
    t=best_t(S)
    if t>=0.5-1e-9:
        print("M",M,"t",t,S); found=True; break
  if found: break
