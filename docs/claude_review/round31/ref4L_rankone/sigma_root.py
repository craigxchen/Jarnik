import numpy as np, sys
from scipy.optimize import root
E=[(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]
inc=[[k for k,e in enumerate(E) if i in e] for i in range(5)]
rng=np.random.default_rng(int(sys.argv[1]))
c=-1.0
def F(s,K):
    r=[]
    for i in range(5):
        r.append(np.prod(s[inc[i]])/K-1)
        r.append(np.sum(s[inc[i]]+c/s[inc[i]]))
    return r
found=[]
for t in range(4000):
    K=np.exp(rng.normal()*2)
    s0=rng.choice([-1,1],size=10)*np.exp(rng.normal(size=10))
    S=root(F,s0,args=(K,),method='hybr')
    if S.success and np.max(np.abs(F(S.x,K)))<1e-11:
        s=S.x; rho=-(s+c/s)*0+(c-s**2)/s
        found.append((np.min(np.diff(np.sort(rho))),K,s))
print(len(found))
nd=[f for f in found if f[0]>1e-6]
print('nondeg',len(nd))
for g,K,s in nd[:10]: print('%.3g K=%.4g'%(g,K),np.round(s,4))
