"""Sequential Monte Carlo estimate of the Haar measure of the gap-admissible set
   Good(WMAX) = { phi in T^s : G_v holds for all v (up to sign) with w(v) <= WMAX },  T = R/(pi/2)Z,
compared with the Poisson prediction exp(-E). Coordinates are added one at a time; at stage j the particles are
extended by a uniform phi_j and filtered by the events whose support has largest index j; the survival fractions
multiply to the estimate (standard SMC with multinomial resampling). Evidence only.
Usage: python3 smc_gaps.py WMAX NPART SEED p1 p2 ..."""
import sys, math, numpy as np
WMAX=float(sys.argv[1]); NP=int(sys.argv[2]); seed=int(sys.argv[3]); ps=[int(x) for x in sys.argv[4:]]
s=len(ps); lw=[math.log(p) for p in ps]; Q=math.pi/2
vs=[]
def dfs(j,v,w,started):
    if j==s:
        if started: vs.append((tuple(v),w))
        return
    m=0
    while w+m*lw[j]<=WMAX+1e-12:
        for sg in ([1] if m==0 else ([1] if not started else [1,-1])):
            v.append(sg*m); dfs(j+1,v,w+m*lw[j],started or m!=0); v.pop()
        m+=1
dfs(0,[],0.0,False)
V=np.array([v for v,w in vs],dtype=float); Wv=np.array([w for v,w in vs])
a=np.arcsin(np.exp(-Wv/2)); b=np.arcsin(np.exp(-Wv/2)/math.sqrt(2)); E=((4/math.pi)*(a+b)).sum()
last=np.array([max(i for i in range(s) if v[i]!=0) for v,w in vs])
rng=np.random.default_rng(seed)
X=np.zeros((NP,s)); logP=0.0
for j in range(s):
    X[:,j]=rng.random(NP)*Q
    idx=np.where(last==j)[0]
    ok=np.ones(NP,bool)
    for c in range(0,len(idx),2000):
        ii=idx[c:c+2000]
        t=np.mod(X@V[ii].T,Q); d0=np.minimum(t,Q-t); d1=np.abs(t-Q/2)
        ok&=~((d0<a[ii][None,:])|(d1<b[ii][None,:])).any(axis=1)
    f=ok.mean()
    if f==0: logP=-math.inf; break
    logP+=math.log(f)
    X=X[ok][rng.integers(0,ok.sum(),NP)]
print("s=%d sum p^-1/2=%.3f #v=%d E=%.3f  log P_good(SMC)=%.3f  -E=%.3f  ratio logP/(-E)=%.3f  primes=%s"%(
    s,sum(p**-0.5 for p in ps),len(vs),E,logP,-E,logP/(-E) if E>0 else float('nan'),ps))
