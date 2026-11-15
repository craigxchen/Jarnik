# evidence-only realisation attempt for each class in a pickle of 6-point rank-5 codes (D=8, c=0)
import numpy as np, pickle, sys
from scipy.optimize import least_squares
cls=list(pickle.load(open(sys.argv[1],'rb')).values()); D=8; L=3; NS=int(sys.argv[2])
rng=np.random.default_rng(7)
for ci,P in enumerate(cls):
    P=np.array(P); d=np.array([(P[i]-P[0])//2 for i in range(1,6)],float)
    def F(v):
        x=v[:D]; y=np.exp(v[D:]); r=x+1j*y
        return np.concatenate([d@np.imag(r**j) for j in range(1,L+1)]+[[x.sum(),y.sum()-D]])
    bestnd=np.inf
    for s in range(NS):
        v0=np.concatenate([rng.normal(size=D)*1.5, rng.normal(size=D)*0.7])
        sol=least_squares(F,v0,xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=3000)
        r=sol.x[:D]+1j*np.exp(sol.x[D:])
        sep=min(abs(r[i]-r[j]) for i in range(D) for j in range(i+1,D))
        if sep>0.05: bestnd=min(bestnd,np.linalg.norm(sol.fun))
    print("class",ci,"best nondegenerate residual (sep>0.05)",bestnd,flush=True)
