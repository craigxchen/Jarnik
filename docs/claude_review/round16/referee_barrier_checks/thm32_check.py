# Numerical test of Thm 3.2 (orbit average = L2 mean) on the sdc7 / sdc11 boxes, incl. the relaxation m < N (N odd).
import sys, itertools, numpy as np, scipy.linalg as sl
sys.path.insert(0,'/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round16/referee_barrier_checks')
from sdc_sym import sdc
def box(q,L):
    P,Pp=sdc(q); pats=[sorted(S) for S in set(P)|set(Pp)]
    K=np.array(list(itertools.product((-1,0,1),repeat=q)),dtype=np.int64)
    ok=np.ones(len(K),bool)
    for S in pats: ok&=np.abs(K[:,S].sum(1))<=L
    return K[ok]
def monos(d,deg):
    out=[]
    for t in range(deg+1):
        for c in itertools.combinations_with_replacement(range(d),t): out.append(c)
    return out
def V(A,m):
    d=A.shape[1]; ms=monos(d,m-1)
    E=np.ones((len(A),len(ms)))
    for i,c in enumerate(ms):
        for x in c: E[:,i]*=A[:,x]
    # annihilators: c with c^T E = 0
    N=sl.null_space(E.T)
    return N
rng=np.random.default_rng(1)
for q,L,ms in [(7,2,[3,4,5,6,7,8]),(11,2,[3,4,5])]:
    A=box(q,L); print(f"sdc{q} L={L} |A|={len(A)}",flush=True)
    # translation g: x->x+1 acting on coords; (g c)_k = c_{g^{-1}k}
    idx={tuple(k):i for i,k in enumerate(A)}
    perm=np.array([idx[tuple(np.roll(k,1))] for k in A])  # k -> shifted k
    for m in ms:
        Nv=V(A,m)
        if Nv.shape[1]==0: print(f"  m={m}: V=0"); continue
        c=Nv@rng.standard_normal(Nv.shape[1])
        a=1; th=np.linspace(0,2*np.pi,4001)[:-1]
        x=np.arange(q)
        def Phi(cc,theta):
            u=np.cos(theta+2*np.pi*a*x/q)
            return (cc*(A@u)**m).sum()
        f=np.array([Phi(c,t) for t in th]); L2=np.sqrt((f**2).mean()); Linf=np.abs(f).max()
        worst=np.inf
        for t0 in th[::97]:
            cs=c.copy(); vals=[]
            for s in range(q):
                vals.append(abs(Phi(cs,t0))); cs=cs[perm]
            # orbit mean of squares vs L2^2
            worst=min(worst,np.mean(np.square(vals))/L2**2)
        print(f"  m={m} dimV={Nv.shape[1]} 2m<q={2*m<q} m<q={m<q}: min_theta0 (orbit-mean |Phi|^2)/||f||_2^2 = {worst:.6f} ; ||f||inf/||f||2={Linf/L2:.3f} sqrt(m+1)={np.sqrt(m+1):.3f}",flush=True)
