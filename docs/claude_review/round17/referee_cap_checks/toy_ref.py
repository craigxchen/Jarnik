# Referee toy check (EVIDENCE, float64): sdc_q box with slab bound L; Psi_A on Fourier circles, N_i splitting (4.1),
# restriction ratio rho_A = mu(A) Psi_A / Psi_T vs permuted circles and uniform random u. Own code.
import numpy as np, itertools, sys, time, scipy.linalg as sla
from math import factorial, comb
def chi(a,q):
    a%=q
    return 0 if a==0 else (1 if pow(a,(q-1)//2,q)==1 else -1)
def box(q,L):
    pats=[[a]+[x for x in range(q) if chi(a-x,q)==1] for a in range(q)]+[[a]+[x for x in range(q) if chi(a-x,q)==-1] for a in range(q)]
    Pm=np.zeros((2*q,q),dtype=np.int64)
    for i,S in enumerate(pats): Pm[i,S]=1
    pts=np.array(list(itertools.product([-1,0,1],repeat=q)),dtype=np.int64)
    ok=np.all(np.abs(pts@Pm.T)<=L,axis=1)
    return pts[ok].astype(float)
def monos(q,maxdeg):
    out=[]
    for e in itertools.product(range(3),repeat=q):
        if sum(e)<=maxdeg: out.append(e)
    return out
def orth_complement_basis(A,m):
    q=A.shape[1]
    ms=[e for e in monos(q,m-1)]
    B=np.ones((A.shape[0],len(ms)))
    for c,e in enumerate(ms):
        for x in range(q):
            if e[x]: B[:,c]*=A[:,x]**e[x]
    Q,R,piv=sla.qr(B,mode='economic',pivoting=True)
    d=np.abs(np.diag(R)); rank=int(np.sum(d>d[0]*1e-10))
    return Q[:,:rank],rank
def psiT(u,m):
    # cube top component norm: (m!)^2 [t^m] prod_x (1 + (2/3)u_x^2 t + (1/4)(2/9)u_x^4 t^2)
    p=np.array([1.0])
    for ux in u: p=np.convolve(p,[1.0,(2/3)*ux**2,(2/9)*ux**4/4])
    return factorial(m)**2*p[m]
def run(q,L,ms,seed=1):
    t0=time.time(); A=box(q,L); U=len(A); muA=U/3**q
    print(f"\n=== sdc{q} L={L}: |A|={U} mu(A)={muA:.4f}",flush=True)
    rng=np.random.default_rng(seed)
    perms=[rng.permutation(q) for _ in range(4)]
    unis=[rng.uniform(-1,1,q) for _ in range(4)]
    for m in ms:
        Qb,rank=orth_complement_basis(A,m)
        dimV=U-rank
        def Psi(u):
            g=(A@u)**m
            r=g-Qb@(Qb.T@g)
            return r@r/U
        def rho(u): return muA*Psi(u)/psiT(u,m)
        circ=[]
        for a in (1,2):
            for th in (0.0,0.3,1.1,2.0):
                u=np.cos(th+2*np.pi*a*np.arange(q)/q)
                circ.append((Psi(u),rho(u)))
        psis=[c[0] for c in circ]; rhos=[c[1] for c in circ]
        prho=[rho(np.cos(0.3+2*np.pi*np.arange(q)/q)[p]) for p in perms]
        urho=[rho(u) for u in unis]
        line=f" m={m} dimV={dimV} Psi_circ rel.spread={(max(psis)-min(psis))/max(psis):.2e} rho_circ=[{min(rhos):.4g},{max(rhos):.4g}] perm=[{min(prho):.4g},{max(prho):.4g}] unif=[{min(urho):.4g},{max(urho):.4g}]"
        # N_i splitting
        if m<q:
            z=np.exp(2j*np.pi*np.arange(q)/q)
            kp=A@z; km=A@np.conj(z)
            N=[]
            for i in range(m+1):
                g=kp**i*km**(m-i)
                r=g-Qb@(Qb.T@g)
                N.append(np.vdot(r,r).real/U)
            pred=sum(comb(m,i)**2*N[i] for i in range(m+1))/4**m
            line+=f" | (4.1): Psi/pred-1={psis[0]/pred-1:.1e}, N_i=N_(m-i) maxrel={max(abs(N[i]-N[m-i])/max(N[i],1e-300) for i in range(m+1)):.1e}, minN/maxN={min(N)/max(N):.3g}"
        print(line,f"[{time.time()-t0:.0f}s]",flush=True)
# validate psiT on the full cube q=5 directly
qq=5; T=np.array(list(itertools.product([-1,0,1],repeat=qq)),dtype=float)
for m in (2,3,4):
    Qb,_=orth_complement_basis(T,m); u=np.random.default_rng(3).uniform(-1,1,qq)
    g=(T@u)**m; r=g-Qb@(Qb.T@g)
    print(f"cube q=5 m={m}: direct Psi_T={r@r/len(T):.10g} formula={psiT(u,m):.10g}")
run(7,2,[3,5,6,7,8])
run(11,2,[3,5,6,7])
