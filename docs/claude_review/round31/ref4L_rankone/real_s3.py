import numpy as np
from scipy.optimize import brentq
rng=np.random.default_rng(0)
def esym(v):
    e=[1,0,0,0,0,0]
    for a in v:
        for k in range(5,0,-1): e[k]+=e[k-1]*a
    return e
cnt=0;pos=0
for t in range(8000):
    y=rng.normal(size=4)
    def g(z):
        v=list(y)+[z]; v=np.array(v)-np.mean(v)  # keep e1=0 by centering
        e=esym(v); return e[2]*e[5]+e[3]*e[4]
    zs=np.linspace(-6,6,61); vals=[g(z) for z in zs]
    for a,b,fa,fb in zip(zs,zs[1:],vals,vals[1:]):
        if fa*fb<0:
            z=brentq(g,a,b); v=np.array(list(y)+[z]); v=v-v.mean(); e=esym(v)
            cnt+=1
            if -e[4]>0 and min(abs(v))>1e-3 and min(np.diff(np.sort(v)))>1e-3: pos+=1
print('real points',cnt,' with -e4(y)>0 :',pos)
