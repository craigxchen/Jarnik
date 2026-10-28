# For contact (mu_inf, mu_0): K = {1..mu_inf-1} U {9-mu_0..8}. Try direct LM and staged (add conditions one at a time, warm start).
import sys, numpy as np
from scipy.optimize import least_squares
from kset import *
exec(open('signpat.py').read().split('print(len(pats))')[0])
def Kof(mi,m0): return list(range(1,mi))+list(range(9-m0,9))
def good(v,K):
    d=diag(v); r=np.linalg.norm(eqs(v,K))
    return r<1e-10 and d['imin']>1e-3 and d['sep']>1e-3 and d['rmax']<1e4 and d['rmin']>1e-4, r, d
def staged_path(mi,m0):
    # order: grow mu_inf first to mi, then mu_0
    path=[]
    for a in range(2,mi+1): path.append((a,0))
    for b in range(1,m0+1): path.append((mi,b))
    return path
if __name__=='__main__':
    mi=int(sys.argv[1]); m0=int(sys.argv[2]); N=int(sys.argv[3]); rng=np.random.default_rng(int(sys.argv[4])); out=open(sys.argv[5],'w')
    gauge=2 if m0==0 else 1
    Kt=Kof(mi,m0); ndir=0; nst=0; sols=[]
    for it in range(N):
        eps=pats[rng.integers(len(pats))]
        # direct
        try:
            x,rn=solve(Kt,eps,rng,max(1,30-gauge-4*len(Kt)))
            ok,r,d=good(unpack(x,eps),Kt)
            if ok: ndir+=1; sols.append(('direct',eps,list(unpack(x,eps))))
        except Exception: pass
        # staged
        try:
            path=staged_path(mi,m0); a0,b0=path[0]; K0=Kof(a0,b0)
            x,rn=solve(K0,eps,rng,max(1,30-2-4*len(K0)))
            ok,r,d=good(unpack(x,eps),K0)
            if not ok: continue
            for (a,b) in path[1:]:
                K=Kof(a,b)
                def F(y):
                    v=unpack(y,eps); return np.concatenate([eqs(v,K),[np.log(abs(v[14]))]])
                sol=least_squares(F,x,method='trf',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=6000)
                ok,r,d=good(unpack(sol.x,eps),K)
                if not ok: break
                x=sol.x
            if ok: nst+=1; sols.append(('staged',eps,list(unpack(x,eps))))
        except Exception: pass
        print(it,'direct',ndir,'staged',nst,file=out,flush=True)
    print('RESULT mu=',(mi,m0),'K=',Kt,'N=',N,'direct_nondeg',ndir,'staged_nondeg',nst,file=out,flush=True)
    for s in sols[:5]: print('SOL',s[0],''.join('+' if e>0 else '-' for e in s[1]),s[2],file=out)
