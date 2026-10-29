# Prop 3.7 test at all vertex classes of Pi_nt (numerical evidence).
# For a vertex v with N v integral: rational curves u_j(s) = lam_j prod_B prod_{i<= N v_B} (s - s_{B,i})^{[j in B]}
# (E_B.C = N v_B), passing through e at 4N parameters.  Expected family dimension 3N+1.  Test: solution found,
# family dimension = 3N+1, and the map (family, s) -> T has rank 4 (dominant)  =>  WTB holds at v (Prop 3.7).
import numpy as np, itertools, sys
from fractions import Fraction
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
Bs=[frozenset(c) for r in range(1,5) for c in itertools.combinations(range(5),r)]
V=np.load('/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round27/split_checks/pi_nt_vertices.npy')
def key(v):
    best=None
    for p in itertools.permutations(range(5)):
        for inv in (False,True):
            w=[]
            for B in Bs:
                PB=frozenset(p[a] for a in B)
                if inv: PB=frozenset(range(5))-PB
                w.append((tuple(sorted(PB)),round(v[Bs.index(B)],6)))
            t=tuple(sorted(w))
            if best is None or t<best: best=t
    return best
classes={}
for v in V:
    k=key(v)
    if k not in classes: classes[k]=v
print("vertex classes under S_5 x iota:",len(classes),flush=True)
def solve_vertex(v,tries=30,maxit=300):
    fr=[Fraction(x).limit_denominator(12) for x in v]
    N=1
    for f in fr: N=N*f.denominator//np.gcd(N,f.denominator)
    mult=[int(f*N) for f in fr]
    roots=[(Bi,B) for Bi,B in enumerate(Bs) for _ in range(mult[Bi])]
    R=len(roots); P=4*N
    inc=np.array([[1 if j in B else 0 for (_,B) in roots] for j in range(5)])
    fixed={0:0.0,1:1.0,2:-1.0}
    free=[i for i in range(R) if i not in fixed]
    nfree=len(free)+4+P
    def unpack(w):
        s=np.zeros(R,dtype=complex)
        for i,val in fixed.items(): s[i]=val
        s[free]=w[:len(free)]; lam=np.concatenate([w[len(free):len(free)+4],[1.0]]); si=w[len(free)+4:]
        return s,lam,si
    def U(s,lam,x):
        return lam*np.prod(np.where(inc==1,(x-s)[None,:],1.0),axis=1)
    def eqs(w):
        s,lam,si=unpack(w)
        out=[]
        for x in si:
            u=U(s,lam,x); out.extend(u[:4]/u[4]-1)
        return np.array(out)
    def jac(w,h=1e-8):
        f0=eqs(w); J=np.zeros((len(f0),nfree),dtype=complex)
        for k in range(nfree):
            e=np.zeros(nfree,dtype=complex); e[k]=h; J[:,k]=(eqs(w+e)-f0)/h
        return J
    for t in range(tries):
        w=rng.normal(size=nfree)+1j*rng.normal(size=nfree)
        ok=False
        with np.errstate(all='ignore'):
            for it in range(maxit):
                f=eqs(w)
                if not np.all(np.isfinite(f)): break
                if np.linalg.norm(f)<1e-11: ok=True; break
                J=jac(w); dw=np.linalg.lstsq(J,-f,rcond=None)[0]
                step=1.0 if np.linalg.norm(dw)<5 else 5/np.linalg.norm(dw)
                w=w+step*dw
        if not ok: continue
        s,lam,si=unpack(w)
        if min(abs(a-b) for a,b in itertools.combinations(list(si),2))<1e-5: continue
        if min(abs(a-b) for a,b in itertools.combinations(list(s),2))<1e-5: continue
        J=jac(w); sv=np.linalg.svd(J,compute_uv=False); rk=int(np.sum(sv>1e-7*sv[0]))
        _,_,Vh=np.linalg.svd(J); K=Vh[rk:].conj().T
        x0=0.31+0.27j
        def y(w,x):
            s,lam,_=unpack(w); u=U(s,lam,x); return u[:4]/u[4]
        h=1e-7; cols=[(y(w+h*K[:,k],x0)-y(w,x0))/h for k in range(K.shape[1])]+[(y(w,x0+h)-y(w,x0))/h]
        D=np.array(cols).T; d=np.linalg.svd(D,compute_uv=False); rD=int(np.sum(d>1e-6*d[0]))
        return dict(N=N,dim=nfree-rk,expected=3*N+1,dominance_rank=rD,smallest_sv=d[min(3,len(d)-1)]/d[0])
    return dict(N=N,dim=None)
if __name__=="__main__":
    only=sys.argv[2] if len(sys.argv)>2 else None
    for i,(k,v) in enumerate(classes.items()):
        vals=sorted(set(np.round(v,6)))
        if only=='onesided' and not set(vals)<= {0.0,1.0}: continue
        if only=='frac' and set(vals)<= {0.0,1.0}: continue
        res=solve_vertex(v)
        print(i,"values",vals,res,flush=True)
