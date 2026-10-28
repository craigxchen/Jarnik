# calibration: same homotopy (4,2)->(4,4) for the doubled k=4 problem (3 columns, 14 places), which has solutions (dbl4.py)
import numpy as np, sys, itertools
from scipy.optimize import least_squares
import dbl4
from fastK import polyjac
MEM=dbl4.MEM; NP=14
KA=[1,2,3,7,8]; KN=[5,6]
def EJ(x,eps,K):
    r=x[:NP]+1j*eps*np.exp(x[NP:]); da=np.ones(NP,complex); db=1j*eps*np.exp(x[NP:])
    F=[];rows=[]
    for j in range(3):
        c,Jc=polyjac(r[MEM[j]]); F.append(c[K].imag)
        blk=np.zeros((len(K),2*NP)); sub=Jc[K]
        blk[:,MEM[j]]=(sub*da[MEM[j]]).imag; blk[:,NP+MEM[j]]=(sub*db[MEM[j]]).imag; rows.append(blk)
    return np.concatenate(F),np.vstack(rows)
def gauge(x,eps):
    i=NP-1; r=x[i]+1j*eps[i]*np.exp(x[2*NP-1]); g=np.log(abs(r)); G=np.zeros(2*NP)
    G[i]=x[i]/abs(r)**2; G[2*NP-1]=(eps[i]*np.exp(x[2*NP-1]))**2/abs(r)**2; return g,G
def diagr(r):
    allv=np.concatenate([r,np.conj(r)]); iu=np.triu_indices(len(allv),1)
    return dict(imin=np.min(np.abs(r.imag)),sep=np.min(np.abs(allv[:,None]-allv[None,:])[iu]),rmax=np.max(np.abs(r)),rmin=np.min(np.abs(r)))
def track(x0,eps,rng,nsl,maxsteps=6000):
    EN0,_=EJ(x0,eps,KN); SL=rng.normal(size=(nsl,2*NP)); s0=SL@x0
    def HJ(z):
        x=z[:2*NP]; tau=z[-1]
        FA,JA=EJ(x,eps,KA); FN,JN=EJ(x,eps,KN); g,G=gauge(x,eps)
        H=np.concatenate([FA,FN-(1-tau)*EN0,[g],SL@x-s0])
        J=np.hstack([np.vstack([JA,JN,G[None,:],SL]),np.concatenate([np.zeros(len(FA)),EN0,[0],np.zeros(nsl)])[:,None]])
        return H,J
    z=np.concatenate([x0,[0.0]]); h=0.01; prev=None
    for step in range(maxsteps):
        H,J=HJ(z); _,_,Vt=np.linalg.svd(J); t=Vt[-1]
        if prev is None:
            if t[-1]<0: t=-t
        elif t@prev<0: t=-t
        zp=z+h*t; ok=False
        for it in range(10):
            H,J=HJ(zp); A=np.vstack([J,t[None,:]]); b=-np.concatenate([H,[t@(zp-z)-h]])
            try: d=np.linalg.solve(A,b)
            except np.linalg.LinAlgError: break
            zp=zp+d
            if np.linalg.norm(d)<1e-11*(1+np.linalg.norm(zp)): ok=True; break
        if ok and np.all(np.isfinite(zp)) and np.linalg.norm(HJ(zp)[0])<1e-9:
            if zp[-1]>=1.0 and z[-1]<1.0:
                w=(1.0-z[-1])/(zp[-1]-z[-1]); zz=np.concatenate([z[:-1]+(zp[:-1]-z[:-1])*w,[1.0]])
                for it in range(30):
                    H,J=HJ(zz); d=np.linalg.lstsq(J[:,:-1],-H,rcond=None)[0]; zz[:-1]+=d
                    if np.linalg.norm(d)<1e-13: break
                return 'reached',zz[:-1],1.0
            prev=t; z=zp; h=min(h*1.4,0.05)
            d=diagr(z[:NP]+1j*eps*np.exp(z[NP:2*NP]))
            if d['sep']<1e-6 or d['imin']<1e-7 or d['rmax']>1e6 or d['rmin']<1e-6: return 'degenerate',z[:-1],z[-1]
            if z[-1]<-0.5: return 'turned_back',z[:-1],z[-1]
        else:
            h/=3
            if h<1e-11: return 'stalled',z[:-1],z[-1]
    return 'maxsteps',z[:-1],z[-1]
if __name__=='__main__':
    N=int(sys.argv[1]); rng=np.random.default_rng(int(sys.argv[2])); out=open(sys.argv[3],'w')
    pats=[e for e in itertools.product([1,-1],repeat=14) if all(sum(e[m] for m in MEM[j])==0 for j in range(3))]
    stats={}
    for it in range(N):
        eps=np.array(pats[rng.integers(len(pats))],float)
        # find a (4,2) solution: 15 eqs (3x5) + gauge, 28 unknowns -> slices 12
        nsl=28-15-1; L=rng.normal(size=(nsl,28))
        def F(x):
            FA,_=EJ(x,eps,KA); g,_=gauge(x,eps); return np.concatenate([FA,[g],0.1*(L@x-L@x0)])
        def Jf(x):
            _,JA=EJ(x,eps,KA); _,G=gauge(x,eps); return np.vstack([JA,G[None,:],0.1*L])
        x0=np.concatenate([rng.normal(size=14)*rng.choice([0.5,1,2]),rng.normal(size=14)*0.7])
        sol=least_squares(F,x0+0.3*rng.normal(size=28),jac=Jf,method='lm',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=2000)
        x=sol.x; r=x[:NP]+1j*eps*np.exp(x[NP:]); d=diagr(r)
        if np.linalg.norm(EJ(x,eps,KA)[0])>1e-10 or d['sep']<1e-3 or d['imin']<1e-3: continue
        # 6 eqs added (3 cols x 2): solution set at tau has dim 28-15-6-1 = 6 ; slices 6
        try: st,xe,tau=track(x,eps,rng,6)
        except Exception: st,xe,tau='error',x,0
        re=xe[:NP]+1j*eps*np.exp(xe[NP:]); de=diagr(re)
        res=np.linalg.norm(np.concatenate([EJ(xe,eps,KA)[0],EJ(xe,eps,KN)[0]]))
        good= st=='reached' and res<1e-8 and de['sep']>1e-3 and de['imin']>1e-3
        st2=st+('-NONDEG' if good else '')
        stats[st2]=stats.get(st2,0)+1
        print(it,st2,f"tau={tau:.4f} res={res:.1e}",' '.join(f"{a}={b:.1e}" for a,b in de.items()),file=out,flush=True)
    print('stats',stats,file=out)
