# complexified Moebius form: g=(g1,g2), gp=(gp1,gp2) in C[s]_{<=4}, g1 gp2 - g2 gp1 = s^4 ; zhat_a in C[s]_{<=1}^2
# z_a = det(zhat_a, g) = zh1 g2 - zh2 g1 ; zp_a = det(zhat_a, gp) ; for each pair: common root rho_ab of z_a,z_b ; rho'_ab of zp_a,zp_b
import numpy as np, itertools, sys
PAIRS=list(itertools.combinations(range(5),2))
NU=20+20+20
def unpack(u):
    g1,g2,h1,h2=u[0:5],u[5:10],u[10:15],u[15:20]   # coefficient arrays high->low (deg 4)
    zh=u[20:40].reshape(5,4)   # zh[a] = (a1,a0,b1,b0): zhat_a = (a1 s + a0, b1 s + b0)
    rho=u[40:50]; rhop=u[50:60]
    return g1,g2,h1,h2,zh,rho,rhop
def zpoly(zh,g1,g2):
    return np.polysub(np.polymul([zh[0],zh[1]],g2),np.polymul([zh[2],zh[3]],g1))
def F(u):
    g1,g2,h1,h2,zh,rho,rhop=unpack(u)
    out=list(np.polysub(np.polymul(g1,h2),np.polymul(g2,h1))-np.array([0,0,0,0,1,0,0,0,0]))
    Z=[zpoly(zh[a],g1,g2) for a in range(5)]; Zp=[zpoly(zh[a],h1,h2) for a in range(5)]
    for k,(a,b) in enumerate(PAIRS):
        out+= [np.polyval(Z[a],rho[k]),np.polyval(Z[b],rho[k]),np.polyval(Zp[a],rhop[k]),np.polyval(Zp[b],rhop[k])]
    return np.array(out,dtype=complex)
def jac(u,h=1e-7):
    f0=F(u); J=np.zeros((len(f0),len(u)),dtype=complex)
    for k in range(len(u)):
        e=np.zeros(len(u),complex); e[k]=h; J[:,k]=(F(u+e)-F(u-e))/(2*h)
    return J
def places(u):
    g1,g2,h1,h2,zh,rho,rhop=unpack(u)
    Z=[zpoly(zh[a],g1,g2) for a in range(5)]
    verts=[]
    for a in range(5):
        rts=np.roots(Z[a]); edges=[rho[k] for k,(p,q) in enumerate(PAIRS) if a in (p,q)]
        rem=[r for r in rts if min(abs(r-e) for e in edges)>1e-6]
        verts.append(rem)
    return rho,rhop,verts
if __name__=='__main__':
    N=int(sys.argv[1]); rng=np.random.default_rng(int(sys.argv[2])); out=open(sys.argv[3],'w'); nd=0; nc=0
    for it in range(N):
        u=rng.normal(size=NU)+1j*rng.normal(size=NU)
        L=rng.normal(size=(11,NU))+1j*rng.normal(size=(11,NU)); b=L@u
        conv=False
        for step in range(200):
            f=np.concatenate([F(u),L@u-b]); 
            J=np.vstack([jac(u),L])
            try: d=np.linalg.solve(J,-f)
            except np.linalg.LinAlgError: break
            lam=1.0; n0=np.linalg.norm(f)
            while lam>1e-4:
                un=u+lam*d
                if np.linalg.norm(np.concatenate([F(un),L@un-b]))<n0*(1-0.25*lam): break
                lam/=2
            u=un
            if np.linalg.norm(np.concatenate([F(u),L@u-b]))<1e-11: conv=True; break
            if not np.all(np.isfinite(u)) or np.max(np.abs(u))>1e7: break
        if not conv: continue
        nc+=1
        rho,rhop,verts=places(u)
        nv=[len(v) for v in verts]
        allp=np.concatenate([rho,rhop]+[np.array(v) for v in verts])
        iu=np.triu_indices(len(allp),1); sep=np.min(np.abs(allp[:,None]-allp[None,:])[iu])
        mn=np.min(np.abs(allp))
        g1,g2,h1,h2,zh,_,_=unpack(u)
        degs=[abs(zpoly(zh[a],g1,g2)[0]) for a in range(5)]
        ok= all(n==1 for n in nv) and sep>1e-4 and mn>1e-4 and np.max(np.abs(allp))<1e4 and min(degs)>1e-6
        nd+=ok
        print(it,'nverts',nv,f"sep={sep:.1e} min|rho|={mn:.1e} max|rho|={np.max(np.abs(allp)):.1e} minlead={min(degs):.1e}",'NONDEG' if ok else '',list(u) if ok else '',file=out,flush=True)
    print('converged',nc,'nondeg',nd,file=out)
