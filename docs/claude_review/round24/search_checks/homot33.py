# path-tracking homotopy from a nondegenerate (4,2) solution to (4,4):
#  H(x,tau) = [E_A(x); E_N(x) - (1-tau) E_N(x0); gauge(x); slice(x)],  E_A: Im e_k, k in {1,2,3,7,8}; E_N: k in {5,6}
import numpy as np, sys, itertools
from scipy.optimize import least_squares
from fastK import Prob, polyjac, MEM, diag, isgood, Kof
exec(open('signpat.py').read().split('print(len(pats))')[0])
import os
KA=[1,2,3,7,8]; KN=[5,6]
if os.environ.get('HSTART')=='33': KA=[1,2,6,7,8]; KN=[3,5]
def EJ(x,eps,K):
    r=x[:15]+1j*eps*np.exp(x[15:]); da=np.ones(15,complex); db=1j*eps*np.exp(x[15:])
    F=[];rows=[]
    for j in range(4):
        c,Jc=polyjac(r[MEM[j]]); F.append(c[K].imag)
        blk=np.zeros((len(K),30)); sub=Jc[K]
        blk[:,MEM[j]]=(sub*da[MEM[j]]).imag; blk[:,15+MEM[j]]=(sub*db[MEM[j]]).imag; rows.append(blk)
    return np.concatenate(F),np.vstack(rows)
def gauge(x,eps):
    r14=x[14]+1j*eps[14]*np.exp(x[29]); g=np.log(abs(r14))
    G=np.zeros(30); G[14]=x[14]/abs(r14)**2; G[29]=(eps[14]*np.exp(x[29]))**2/abs(r14)**2
    return g,G
def Hfun(x,tau,eps,EN0,sl,s0):
    FA,JA=EJ(x,eps,KA); FN,JN=EJ(x,eps,KN); g,G=gauge(x,eps)
    H=np.concatenate([FA,FN-(1-tau)*EN0,[g],[sl@x-s0]])
    Jx=np.vstack([JA,JN,G[None,:],sl[None,:]]); Jt=np.concatenate([np.zeros(len(FA)),EN0,[0,0]])
    return H,Jx,Jt
def track(x0,eps,rng,out,maxsteps=6000):
    EN0,_=EJ(x0,eps,KN); sl=rng.normal(size=30); s0=sl@x0
    z=np.concatenate([x0,[0.0]]); h=0.01; prev=None
    def HJ(z):
        H,Jx,Jt=Hfun(z[:30],z[30],eps,EN0,sl,s0); return H,np.hstack([Jx,Jt[:,None]])
    for step in range(maxsteps):
        H,J=HJ(z)
        _,_,Vt=np.linalg.svd(J); t=Vt[-1]
        if prev is None:
            if t[30]<0: t=-t
        elif t@prev<0: t=-t
        zp=z+h*t; ok=False
        for it in range(10):
            H,J=HJ(zp)
            A=np.vstack([J,t[None,:]]); b=-np.concatenate([H,[t@(zp-z)-h]])
            try: d=np.linalg.solve(A,b)
            except np.linalg.LinAlgError: break
            zp=zp+d
            if np.linalg.norm(d)<1e-11*(1+np.linalg.norm(zp)): ok=True; break
        if ok and np.all(np.isfinite(zp)) and np.linalg.norm(HJ(zp)[0])<1e-9:
            if zp[30]>=1.0 and z[30]<1.0:
                x=z[:30]+(zp[:30]-z[:30])*(1.0-z[30])/(zp[30]-z[30])
                for it in range(20):
                    Hh,Jx,Jt=Hfun(x,1.0,eps,EN0,sl,s0); d=np.linalg.lstsq(Jx,-Hh,rcond=None)[0]; x=x+d
                    if np.linalg.norm(d)<1e-13: break
                return 'reached',x,1.0
            prev=t; z=zp; h=min(h*1.4,0.05)
            r=z[:15]+1j*eps*np.exp(z[15:30]); d=diag(r)
            if d['sep']<1e-6 or d['imin']<1e-7 or d['rmax']>1e6 or d['rmin']<1e-6: return 'degenerate',z[:30],z[30]
            if z[30]<-0.5: return 'turned_back',z[:30],z[30]
        else:
            h/=3
            if h<1e-11: return 'stalled',z[:30],z[30]
    return 'maxsteps',z[:30],z[30]
if __name__=='__main__':
    N=int(sys.argv[1]); rng=np.random.default_rng(int(sys.argv[2])); out=open(sys.argv[3],'w')
    stats={}
    for it in range(N):
        eps=np.array(pats[rng.integers(len(pats))],float)
        # stage to (4,2) by staged LM
        x=None
        for K in (([1,2,6,7,8],) if os.environ.get('HSTART')=='33' else ([1,2,3],[1,2,3,8],[1,2,3,7,8])):
            ns=max(1,29-4*len(K)) if x is None else 0
            P=Prob(K,eps,max(ns,1),rng)
            if x is None:
                x0=np.concatenate([rng.normal(size=15)*rng.choice([0.5,1,2]), rng.normal(size=15)*0.7])
                sol=least_squares(P.F,x0,jac=P.J,method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=2000)
            else:
                f=lambda y: P.F(y)[:4*len(K)+1]; jf=lambda y: P.J(y)[:4*len(K)+1]
                sol=least_squares(f,x,jac=jf,method='trf',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=3000)
            x=sol.x; r=P.rho(x)
            if not isgood(r,P.eqnorm(x)): x=None; break
        if x is None: continue
        try: status,xe,tau=track(x,eps,rng,out)
        except Exception as e: status,xe,tau='error',x,0.0
        stats[status]=stats.get(status,0)+1
        r=xe[:15]+1j*eps*np.exp(xe[15:]); d=diag(r)
        res=np.linalg.norm(np.concatenate([EJ(xe,eps,KA)[0],EJ(xe,eps,KN)[0]]))
        FULL=frozenset(range(1,6)); from fastK import SUBS as SB
        PL=list(SB)+[FULL-T for T in SB]; allv=np.concatenate([r,np.conj(r)]); lab=lambda A:''.join(str(a) for a in sorted(A))
        coinc=[(lab(PL[a]),lab(PL[b])) for a in range(30) for b in range(a+1,30) if abs(allv[a]-allv[b])<1e-3] if tau>0.99 else []
        realp=[lab(PL[a]) for a in range(15) if abs(r[a].imag)<1e-3] if tau>0.99 else []
        print(it,status,'coinc',coinc,'real',realp,file=out)
        print(it,status,f"tau={tau:.4f} res44={res:.1e}",' '.join(f"{a}={b:.1e}" for a,b in d.items()), list(r) if tau>0.99 else '',file=out,flush=True)
    print('stats',stats,file=out)
