# Exact: (S_min) fails for every n. u=(1, X_1+i,...,X_{n-1}+i, (X_1+i)(X_2+i)), X_x=T+x.
import sys; sys.path.insert(0,'.')
from gauss import *
from itertools import combinations
def check(name,u,rank):
    assert nonreal_pairs(u)
    n=len(u)-1; hl,G=h_line(u); s,a=minanchor(u)
    print(f"{name}: n={n} rank<={rank} logT={log(T):.2f} h(ell)={hl:.3f} (2logT={2*log(T):.3f}) G={G} "
          f"minanchor={s:.3f} ratio/n={s/hl/n:.4f}  S-(n/3)h(ell)={s-n/3*hl:.2f}",flush=True)
for T in [10**10,10**30,10**60]:
  for n in [2,3,4,6,10,20,40]:
    v=[(T+x,1) for x in range(1,n)]
    u=[(1,0)]+v+[mul(v[0],v[1])] if n>=3 else [(1,0),(T+1,1),pw((T+1,1),2)]
    check(f"one-relation T=1e{len(str(T))-1}",u,n-1)
T=10**30
for k in [4,6,8,10]:
    v=[(T+x,1) for x in range(1,k+1)]
    u=[(1,0)]+v+[mul(v[a],v[b]) for a,b in combinations(range(k),2)]
    check(f"pair-products k={k}",u,k)
