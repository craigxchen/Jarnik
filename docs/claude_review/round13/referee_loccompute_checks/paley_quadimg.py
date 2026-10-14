import numpy as np
from paley_base import build
from rlib import J_exact
A=build(2)
for S in [(0,1,5),(0,2,6)]:
    U=np.array([[1]*11]+[[1 if j==i else 0 for j in range(11)] for i in S])
    img=sorted(set(map(tuple,(A@U.T).tolist())))
    J=J_exact(img,9)
    H=[]; 
    from math import comb
    for D in range(10): H.append(comb(D+3,3)-len(J[D]))
    print("sigma+coords",S,"|img|",len(img),"dim J_D",[len(J[D]) for D in range(10)],"H",H,"cum",np.cumsum(H).tolist())
    full=[(s,a,b,c) for s in range(-2,3) for a in (-1,0,1) for b in (-1,0,1) for c in (-1,0,1)]
    print("  missing from product:",sorted(set(full)-set(img)))
