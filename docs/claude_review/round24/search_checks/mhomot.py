# homotopy in the Moebius/M form: M real 2x2 poly matrix deg<=4 with det M = s^4 ; zhat_a real linear 2-vectors;
# z_a = (1,i) M zhat_a ; rho_ab = chosen coincidence point of zhat_a, zhat_b ; equations E_ab = z_a(rho_ab) (complex) = 0
import numpy as np, itertools, sys
PAIRS=list(itertools.combinations(range(5),2))
S4=np.array([0,0,0,0,1,0,0,0,0],float)
def unpackM(u): return u[:20].reshape(2,2,5)   # M[i,j] coeffs high->low deg 4
def unpackZ(u): return u[20:40].reshape(5,2,2) # zhat[a][i] = (coef s, coef 1)
def detM(M): return np.polysub(np.polymul(M[0,0],M[1,1]),np.polymul(M[0,1],M[1,0]))
def zpoly(M,zh):
    # (1,i) M zh : row vector (M00 + i M10, M01 + i M11) dot (zh0, zh1)
    l0=M[0,0]+1j*M[1,0]; l1=M[0,1]+1j*M[1,1]
    return np.polyadd(np.polymul(l0,zh[0]),np.polymul(l1,zh[1]))
def coinc(zha,zhb,sg):
    q=np.polysub(np.polymul(zha[0],zhb[1]),np.polymul(zha[1],zhb[0]))  # real quadratic
    r=np.roots(q)
    if len(r)<2 or abs(r[0].imag)<1e-14: return None
    return r[0] if (r[0].imag>0)==(sg>0) else r[1]
def E(u,signs):
    M=unpackM(u); Z=unpackZ(u); out=[]
    for (a,b),sg in zip(PAIRS,signs):
        rho=coinc(Z[a],Z[b],sg)
        if rho is None: return None
        v=np.polyval(zpoly(M,Z[a]),rho); out+=[v.real,v.imag]
    return np.array(out)
def Hfun(u,tau,signs,E0,SL,s0):
    e=E(u,signs)
    if e is None: return None
    d=detM(unpackM(u)); d=np.concatenate([np.zeros(9-len(d)),d])-S4
    return np.concatenate([d,e-(1-tau)*E0,SL@u-s0])
def jac(u,tau,signs,E0,SL,s0,h=1e-7):
    f0=Hfun(u,tau,signs,E0,SL,s0); J=np.zeros((len(f0),len(u)))
    for k in range(len(u)):
        e=np.zeros(len(u)); e[k]=h
        fp=Hfun(u+e,tau,signs,E0,SL,s0); fm=Hfun(u-e,tau,signs,E0,SL,s0)
        if fp is None or fm is None: return None,None
        J[:,k]=(fp-fm)/(2*h)
    Jt=np.concatenate([np.zeros(9),E0,np.zeros(len(s0))])
    return J,Jt
def random_start(rng):
    while True:
        a=rng.normal(size=5); d=rng.normal(size=5)
        oct_=np.polysub(np.polymul(a,d),S4)   # = b*c
        r=np.roots(oct_)
        # split roots into two conjugation-closed groups of 4
        rr=sorted(r,key=lambda z:(abs(z.imag)<1e-9, z.real))
        cplx=[z for z in r if z.imag>1e-9]; real=sorted([z.real for z in r if abs(z.imag)<=1e-9])
        groups=[[],[]]
        items=[[z,np.conj(z)] for z in cplx]+[[x] for x in real]
        rng.shuffle(items)
        for it in items:
            g=0 if len(groups[0])+len(it)<=4 and (len(groups[1])>=4 or rng.random()<0.5) else 1
            if len(groups[g])+len(it)>4: g=1-g
            if len(groups[g])+len(it)>4: break
            groups[g]+=it
        if len(groups[0])!=4 or len(groups[1])!=4: continue
        lc=oct_[0]
        b=np.real(np.poly(groups[0]))*lc; c=np.real(np.poly(groups[1]))
        M=np.array([[a,b],[c,d]])
        if np.max(np.abs(np.polysub(detM(M),S4)[-9:]))>1e-8: continue
        # zhat: rotations about random points (pairwise elliptic)
        from mob import sample_elliptic
        Ms=sample_elliptic(rng)
        Z=np.array([[[Mm[0,0],Mm[0,1]],[Mm[1,0],Mm[1,1]]] for Mm in Ms])
        return np.concatenate([M.reshape(-1),Z.reshape(-1)])
def places(u,signs):
    M=unpackM(u); Z=unpackZ(u); rhos=[coinc(Z[a],Z[b],sg) for (a,b),sg in zip(PAIRS,signs)]
    verts=[]; degs=[]
    for a in range(5):
        zp=zpoly(M,Z[a]); degs.append(abs(zp[0])); rts=np.roots(zp)
        edges=[rhos[k] for k,(p,q) in enumerate(PAIRS) if a in (p,q)]
        verts.append([r for r in rts if min(abs(r-e) for e in edges)>1e-6])
    return rhos,verts,degs
def track(u0,signs,rng,maxsteps=700):
    E0=E(u0,signs); SL=rng.normal(size=(11,40)); s0=SL@u0
    z=np.concatenate([u0,[0.0]]); h=0.01; prev=None
    def HJ(z):
        H=Hfun(z[:40],z[40],signs,E0,SL,s0)
        if H is None: return None,None
        J,Jt=jac(z[:40],z[40],signs,E0,SL,s0)
        if J is None: return None,None
        return H,np.hstack([J,Jt[:,None]])
    for step in range(maxsteps):
        H,J=HJ(z)
        if H is None: return 'undefined',z[:40],z[40]
        _,_,Vt=np.linalg.svd(J); t=Vt[-1]
        if prev is None:
            if t[-1]<0: t=-t
        elif t@prev<0: t=-t
        zp=z+h*t; ok=False
        for it in range(8):
            H,J=HJ(zp)
            if H is None: break
            A=np.vstack([J,t[None,:]]); b=-np.concatenate([H,[t@(zp-z)-h]])
            try: d=np.linalg.solve(A,b)
            except np.linalg.LinAlgError: break
            zp=zp+d
            if np.linalg.norm(d)<1e-10*(1+np.linalg.norm(zp)): ok=True; break
        if ok and np.linalg.norm(HJ(zp)[0])<1e-8:
            if zp[-1]>=1.0 and z[-1]<1.0:
                w=(1-z[-1])/(zp[-1]-z[-1]); zz=np.concatenate([z[:-1]+w*(zp[:-1]-z[:-1]),[1.0]])
                for it in range(30):
                    H,J=HJ(zz)
                    if H is None: break
                    dd=np.linalg.lstsq(J[:,:-1],-H,rcond=None)[0]; zz[:-1]+=dd
                    if np.linalg.norm(dd)<1e-12: break
                return 'reached',zz[:-1],1.0
            prev=t; z=zp; h=min(h*1.4,0.05)
            if np.max(np.abs(z))>1e5: return 'diverged',z[:-1],z[-1]
            if z[-1]<-0.5: return 'turned_back',z[:-1],z[-1]
        else:
            h/=3
            if h<1e-10: return 'stalled',z[:-1],z[-1]
    return 'maxsteps',z[:-1],z[-1]
if __name__=='__main__':
    N=int(sys.argv[1]); rng=np.random.default_rng(int(sys.argv[2])); out=open(sys.argv[3],'w'); stats={}
    for it in range(N):
        u0=random_start(rng); signs=rng.choice([1,-1],size=10)
        if E(u0,signs) is None: continue
        try: st,u,tau=track(u0,signs,rng)
        except Exception as e: st,u,tau='error:'+type(e).__name__,u0,0
        info=''
        if st=='reached':
            rhos,verts,degs=places(u,signs)
            nv=[len(v) for v in verts]
            if all(n==1 for n in nv) and min(degs)>1e-8:
                allp=np.array(list(rhos)+[v[0] for v in verts]); allv=np.concatenate([allp,np.conj(allp)])
                iu=np.triu_indices(30,1); sep=np.min(np.abs(allv[:,None]-allv[None,:])[iu]); imin=np.min(np.abs(allp.imag))
                good=sep>1e-4 and imin>1e-4 and np.min(np.abs(allp))>1e-5 and np.max(np.abs(allp))<1e5
                info=f"sep={sep:.1e} imin={imin:.1e} minlead={min(degs):.1e}"
                if good: st='reached-NONDEG'; info+=' '+repr((list(u),list(signs)))
            else: info=f"nverts={nv} minlead={min(degs):.1e}"
        stats[st]=stats.get(st,0)+1
        print(it,st,f"tau={tau:.4f}",info,file=out,flush=True)
    print('stats',stats,file=out)
