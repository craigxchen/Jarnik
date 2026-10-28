import sys, numpy as np, time
from fastK import *
exec(open('signpat.py').read().split('print(len(pats))')[0])
mi,m0,N,seed=map(int,sys.argv[1:5])
nfev=int(sys.argv[6]) if len(sys.argv)>6 else 3000; out=open(sys.argv[5],'w'); rng=np.random.default_rng(seed)
K=Kof(mi,m0); nd=0; nconv=0; t0=time.time(); sols=[]
for it in range(N):
    eps=pats[rng.integers(len(pats))] if (mi,m0) in [(4,4),(8,0)] else tuple(rng.choice([1,-1],size=15))
    try: x,res,P=run_direct(K,eps,rng,2 if m0==0 else 1,nfev=nfev)
    except Exception: continue
    r=P.rho(x)
    if res<1e-10: nconv+=1
    if isgood(r,res): nd+=1; sols.append(list(r))
print('mu',(mi,m0),'K',K,'N',N,'converged',nconv,'nondeg',nd,'time',round(time.time()-t0,1),file=out,flush=True)
for s in sols[:3]: print('SOL',s,file=out)
