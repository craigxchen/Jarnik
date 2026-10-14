import numpy as np
from paley_base import build, COLS
from rlib import J_exact, monos
A=build(2)
y0=COLS[0]; s=np.ones(11,dtype=np.int64); s[2]-=1
img=sorted(set(map(tuple,np.stack([A@y0,A@s],axis=1).tolist())))
print("image",len(img)); print(img)
J=J_exact(img,10)
print("dim J(img)_D D<=10",[len(J[D]) for D in range(11)])
# rho: least D with sum_{i<=D} H(i) = |img|; H(D)=D+1-dimJ_D
H=[D+1-len(J[D]) for D in range(11)]; cum=np.cumsum(H); print("H",H,"cum",cum.tolist())
