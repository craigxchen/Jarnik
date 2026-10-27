# Evidence: simulated annealing over ternary f (F=sum f arbitrary in {-1,0,1} enforced by rejection) minimising sigma.
# Records best sigma found per l1-norm bucket. Exact recomputation of every reported value.
import numpy as np, sys, time
from hlib import *
def run(q,iters,seed,T0=3.0,init=None):
    rng=np.random.default_rng(seed); ch=chi_vec(q)
    G=np.array([np.roll(ch,x) for x in range(q)])  # G[x,a]=chi(a-x)
    if init is None:
        f=np.zeros(q,dtype=np.int64); k=rng.integers(2,q//2)
        k=2*(k//2)+2; idx=rng.choice(q,size=k,replace=False); f[idx[:k//2]]=1; f[idx[k//2:]]=-1
    else: f=init.copy()
    F=int(f.sum()); c=f@G; 
    def Nval(f,F,c): return int(np.maximum(np.abs(F+f),np.abs(c)).sum())
    def ok(f,F):
        s=int(np.abs(f).sum())
        if abs(F)>1: return False
        if s==0: return False
        if F==0 and s==2: return False
        if abs(F)==1 and s==1: return False
        return True
    cur=Nval(f,F,c) if ok(f,F) else 10**9
    best={}
    for it in range(iters):
        T=T0*(1-it/iters)+0.05
        x=rng.integers(q); nv=rng.choice([v for v in (-1,0,1) if v!=f[x]])
        d=nv-f[x]; F2=F+d
        f[x]=nv; c2=c+d*G[x]
        if ok(f,F2):
            v=Nval(f,F2,c2)
            if v<=cur or rng.random()<np.exp((cur-v)/T):
                cur=v; F=F2; c=c2
                s=int(np.abs(f).sum())+abs(F)
                b=s
                if b not in best or v<best[b][0]: best[b]=(v,f.copy())
                continue
        f[x]-=d
    return best
if __name__=="__main__":
    q=int(sys.argv[1]); iters=int(sys.argv[2]); seeds=int(sys.argv[3])
    allbest={}
    t0=time.time()
    for sd in range(seeds):
        b=run(q,iters,sd)
        for k,(v,f) in b.items():
            if k not in allbest or v<allbest[k][0]: allbest[k]=(v,f)
        gm=min(v for v,_ in allbest.values())
        print(f"seed {sd} done t={time.time()-t0:.0f}s global min sigma={gm-q} sigma/r={(gm-q)/(2*q):.4f}",flush=True)
    # buckets
    ks=sorted(allbest)
    print("||lambda||_1 bucket : min sigma/r found (exact recheck)")
    ch=chi_vec(q)
    for k in ks:
        v,f=allbest[k]
        s=sigma_formula(f,q,ch); assert s==v-q
    bk={}
    for k in ks:
        B=min(k//max(1,q//20),40)
        v=allbest[k][0]-q
        if B not in bk or v<bk[B][0]: bk[B]=(v,k)
    for B in sorted(bk): print(f"  l1 in bucket {B} (rep l1={bk[B][1]}): min sigma/r={bk[B][0]/(2*q):.4f}")
