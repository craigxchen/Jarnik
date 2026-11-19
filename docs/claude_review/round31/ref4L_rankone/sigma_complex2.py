import numpy as np, sys
from scipy.optimize import root
E=[(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]
inc=[[k for k,e in enumerate(E) if i in e] for i in range(5)]
rng=np.random.default_rng(int(sys.argv[1])); c=-1.0
def Fc(z,K):
    s=z[:10]+1j*z[10:]; r=[]
    for i in range(5):
        r.append(np.prod(s[inc[i]])/K-1); r.append(np.sum(s[inc[i]]+c/s[inc[i]]))
    r=np.array(r); return np.concatenate([r.real,r.imag])
conv=0; nd=0; classes={}
for t in range(1500):
    K=complex(rng.normal(),rng.normal())*rng.choice([0.3,1,3])
    S=root(Fc,rng.normal(size=20)*rng.choice([0.5,1,2]),args=(K,),method='hybr')
    if np.max(np.abs(Fc(S.x,K)))<1e-11:
        conv+=1; s=S.x[:10]+1j*S.x[10:]
        rho=(c-s**2)/s
        nd_rho=len(set(np.round(rho,6)))
        classes[nd_rho]=classes.get(nd_rho,0)+1
print('converged',conv,' number of distinct rho values among solutions:',classes)
