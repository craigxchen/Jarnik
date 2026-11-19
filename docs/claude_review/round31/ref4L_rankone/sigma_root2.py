import numpy as np, sys, itertools
from scipy.optimize import root
E=[(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]
inc=[[k for k,e in enumerate(E) if i in e] for i in range(5)]
# Eulerian sign patterns (all vertex degrees even), up to S5 we just take all
pats=[]
for mask in range(1<<10):
    sg=np.array([(-1 if mask>>k&1 else 1) for k in range(10)])
    if all(sum(1 for k in inc[i] if sg[k]<0)%2==0 for i in range(5)): pats.append(sg)
print('patterns',len(pats))
rng=np.random.default_rng(int(sys.argv[1]))
def F(th,sg,kap):
    s=sg*np.exp(th); r=[]
    for i in range(5):
        r.append(np.sum(th[inc[i]])-kap)
        r.append(np.sum(s[inc[i]]-1/s[inc[i]]))
    return r
best=[];cnt=0
for sg in pats:
    for t in range(60):
        kap=rng.normal()*3
        S=root(F,rng.normal(size=10)*1.5,args=(sg,kap),method='lm')
        if np.max(np.abs(F(S.x,sg,kap)))<1e-10:
            cnt+=1
            s=sg*np.exp(S.x); rho=(-1-s**2)/s
            best.append((np.min(np.diff(np.sort(rho))),kap,s))
print('solutions',cnt,' nondeg',sum(1 for b in best if b[0]>1e-5))
for b in sorted(best,key=lambda b:-b[0])[:8]: print('%.3g kap=%.3f'%(b[0],b[1]),np.round(b[2],4))
