import numpy as np
from scipy.optimize import least_squares
P=np.array([[1,1,1,1,1,-1,-1,-1],[1,1,1,-1,-1,1,1,-1],[1,1,-1,1,-1,1,-1,1],[1,-1,1,1,-1,-1,1,1],[1,-1,-1,-1,1,1,1,1],[-1,1,1,-1,1,1,-1,1]])
D=8;L=3
d=np.array([(P[i]-P[0])//2 for i in range(1,6)],float)
rng=np.random.default_rng(1)
def F(v):
    x=v[:D]; y=np.exp(v[D:])  # Im>0
    r=x+1j*y
    eq=[d@np.imag(r**j) for j in range(1,L+1)]
    eq=np.concatenate(eq)
    g=[x.sum(), y.sum()-D]
    return np.concatenate([eq,g])
out=[]
for s in range(400):
    v0=np.concatenate([rng.normal(size=D)*1.5, rng.normal(size=D)*0.7])
    sol=least_squares(F,v0,xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
    x=sol.x[:D];y=np.exp(sol.x[D:]); r=x+1j*y
    sep=min(abs(r[i]-r[j]) for i in range(D) for j in range(i+1,D))
    out.append((np.linalg.norm(sol.fun),sep,y.min()))
out=np.array(out)
good=out[(out[:,0]<1e-10)]
print("starts",len(out),"res<1e-10:",len(good))
if len(good): print("max sep among good",good[:,1].max(),"min Im",good[:,2].max())
o=out[np.argsort(out[:,0])][:8]; print(o)
# best with sep>0.05
nd=out[out[:,1]>0.05]; print("best residual with sep>0.05:", nd[:,0].min() if len(nd) else None)
