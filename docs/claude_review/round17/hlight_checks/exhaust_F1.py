# Exhaustive exact check of the |F|=1 bounds for small supports (lambda_inf=-1; negation covers +1):
# n=||lambda||_1=4: f = d_x + d_y - d_z  (AGL: x=0,y=1)  -> claim sigma >= (q-1)/2
# n=6: f = d_0 + d_1 + d_y - d_z - d_w  (AGL: two +1's at 0,1) -> claim sigma >= (q+1)/2
import numpy as np, itertools
from hlib import *
for q in [7,11,19,23,31,43,47,59,67,71,79,83]:
    ch=chi_vec(q); m4=10**9; m6=10**9
    for z in range(2,q):
        f=np.zeros(q,dtype=np.int64); f[0]=1; f[1]=1; f[z]=-1
        m4=min(m4,sigma_formula(f,q,ch))
    if q<=47:
        for y in range(2,q):
            for z,w in itertools.combinations([t for t in range(2,q) if t!=y],2):
                f=np.zeros(q,dtype=np.int64); f[0]=1; f[1]=1; f[y]=1; f[z]=-1; f[w]=-1
                m6=min(m6,sigma_formula(f,q,ch))
    print(f"q={q} |F|=1 n=4: min sigma={m4} bound (q-1)/2={(q-1)//2} ok={m4>=(q-1)//2};  n=6: min sigma={m6 if q<=47 else 'skip'} bound (q+1)/2={(q+1)//2}",flush=True)
