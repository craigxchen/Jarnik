# Referee pass 2 (EVIDENCE, float64): toy check of Prop 4.1 (exact L^2 splitting, constancy on circles for m<q)
# and of the toy-table claim "for m<q circles have the largest restriction ratio". Own code.
import itertools, sys, numpy as np, scipy.linalg as sl
from math import comb, factorial
out=[]
def P(*a):
    s=" ".join(map(str,a)); print(s); out.append(s); sys.stdout.flush()
def chi(a,q):
    a%=q
    return 0 if a==0 else (1 if pow(a,(q-1)//2,q)==1 else -1)
def box(q,L):
    pats=[[a]+[x for x in range(q) if chi(a-x,q)==s] for a in range(q) for s in (1,-1)]
    M=np.zeros((2*q,q),dtype=np.int64)
    for i,p in enumerate(pats): M[i,p]=1
    T=np.array(list(itertools.product((-1,0,1),repeat=q)),dtype=np.int64)
    keep=np.all(np.abs(T@M.T)<=L,axis=1)
    return T[keep]
def basis(A,m):
    q=A.shape[1]; cols=[]
    for beta in itertools.product(range(3),repeat=q):
        if sum(beta)<m: cols.append(np.prod(A.astype(float)**np.array(beta),axis=1))
    B=np.array(cols).T
    Qm,R,piv=sl.qr(B,mode='economic',pivoting=True)
    rk=int(np.sum(np.abs(np.diag(R))>1e-9*abs(R[0,0])))
    return Qm[:,:rk]
def psiT(u,m):
    poly=np.zeros(m+1); poly[0]=1
    for ux in u:
        new=poly.copy(); new[1:]+= (2/3)*ux**2*poly[:-1]; new[2:]+= (1/18)*ux**4*poly[:-2]; poly=new
    return factorial(m)**2*poly[m]
rng=np.random.default_rng(7171)
for q,L,ms in [(7,2,[3,4,5,6,7]),(11,2,[3,5,6,7])]:
    A=box(q,L); U=len(A); mu=U/3**q
    P(f"\n### sdc{q}, L={L}: |A|={U}")
    zeta=np.exp(2j*np.pi/q)
    kh1=A@zeta**np.arange(q); khm=np.conj(kh1)
    for m in ms:
        Qb=basis(A,m)
        top=lambda f: f-Qb@(Qb.conj().T@f)
        def Psi(u):
            f=(A@u).astype(float)**m; t=top(f); return np.vdot(t,t).real/U
        N=[np.vdot(top(kh1**i*khm**(m-i)),top(kh1**i*khm**(m-i))).real/U for i in range(m+1)]
        rhs=sum(comb(m,i)**2*N[i] for i in range(m+1))/4**m
        circ=[]
        for a in (1,2):
            for th in (0,0.3,1.1,2.0):
                u=np.cos(th+2*np.pi*a*np.arange(q)/q); circ.append((Psi(u),mu*Psi(u)/psiT(u,m)))
        cp=np.array([c[0] for c in circ])
        perm=[]; unif=[]
        for _ in range(6):
            u=np.cos(rng.uniform(0,6.3)+2*np.pi*rng.permutation(q)/q); perm.append(mu*Psi(u)/psiT(u,m))
            u=rng.uniform(-1,1,q); unif.append(mu*Psi(u)/psiT(u,m))
        rc=[c[1] for c in circ]
        P(f"m={m} dimV={U-Qb.shape[1]}: Psi_A circle spread rel={np.ptp(cp)/cp.mean():.1e}, (4.1) rel err={abs(rhs-cp[0])/cp[0]:.1e}, N_i=N_(m-i): {np.allclose(N,N[::-1],rtol=1e-8)}; "
          f"rho circ {min(rc):.4g}-{max(rc):.4g} | perm {min(perm):.4g}-{max(perm):.4g} | unif {min(unif):.4g}-{max(unif):.4g}")
open(__file__.replace('.py','_out.txt'),'w').write("\n".join(out)+"\n")
