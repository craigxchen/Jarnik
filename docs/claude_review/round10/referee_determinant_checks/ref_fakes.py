# Referee: small-scale test of Corollary 4.3 on the two fake families.
#  (1) exact rank of the exponent profile (twins / S'_a - S_a) -> full rank M-1 on the slice, so the
#      low-rank automatic membership of Prop 5.2 does not apply to the fakes;
#  (2) slab-box support A (median slab cut, as in the author's evidence scripts) built from the fake's
#      OWN columns beta_p, exact integer basis of V(A,m);
#  (3) for 200 random fake positions delta ~ U[0,Delta]^M (the construction's distribution), the relation
#      (4.2) fails: some c has |sum_k c_k e^{i<k,delta>}| > 0, certified by the leading Taylor term.
# Sizes are far below the regime n-1 >= 300 ln 4r; this only illustrates Prop 4.2 / Cor 4.3.
import numpy as np, itertools, math, sys
import sympy as sp, mpmath as mp
mp.mp.dps=50
rng=np.random.default_rng(5)
def legendre(a,q):
    a%=q
    if a==0: return 0
    return 1 if pow(a,(q-1)//2,q)==1 else -1
def paley_I(q):
    # rows/cols indexed by {inf} u F_q; H(inf,.)=1, H(i,inf)=-1, H(i,j)=delta_ij+chi(j-i)
    idx=['inf']+list(range(q)); M=q+1
    H=np.zeros((M,M),dtype=int)
    for a,x in enumerate(idx):
        for b,y in enumerate(idx):
            if x=='inf': H[a,b]=1
            elif y=='inf': H[a,b]=-1
            else: H[a,b]=(1 if x==y else 0)+legendre(y-x,q)
    assert (H@H.T==M*np.eye(M,dtype=int)).all()
    return H
def skew_double_core(q):
    K=paley_I(q)  # K + K^T = 2I after normalisation?  build skew form K=I+S
    M=q+1
    # make K skew-Hadamard: Paley-I matrix is I+S with S skew; normalise K(x,0)=1 via D K D
    Kt=K.copy()
    if not ((Kt+Kt.T)==2*np.eye(M,dtype=int)).all():
        Kt=Kt.T
    assert ((Kt+Kt.T)==2*np.eye(M,dtype=int)).all()
    Dg=np.diag(Kt[:,0]); Kt=Dg@Kt@Dg
    assert (Kt[:,0]==1).all() and ((Kt+Kt.T)==2*np.eye(M,dtype=int)).all()
    cols=[]
    for a in range(1,M):
        S=[x for x in range(M) if Kt[x,a]==-1]
        Sp=set(S)^{0,a}
        cols.append([1 if x in S else 0 for x in range(M)])
        cols.append([1 if x in Sp else 0 for x in range(M)])
    return np.array(cols).T   # M x r matrix of a_{p,x}
def round5_profile(q,b=8):
    H=paley_I(q); M=q+1
    Hn=H.copy(); Hn[1:,:]*=-1      # multiply rows i in F_q by -1
    # identity pairing pi(i)=i: label of row i (index 1+i) is column 1+i; x_inf gets label a*=1
    lab={0:1}
    for i in range(q): lab[1+i]=1+i
    cols=[]
    for a in range(1,M):
        for _ in range(b): cols.append(Hn[:,a].copy())
    for x in range(M):
        c=Hn[:,lab[x]].copy(); c[x]*=-1; cols.append(c)
    S=np.array(cols).T
    return (1+S)//2
out=[]
def pr(s): print(s,flush=True); out.append(s)
for name,Aexp in [("skew double core q=7",skew_double_core(7)),("round-5 Paley-I b=8 q=7",round5_profile(7))]:
    M,r=Aexp.shape
    diffs=sp.Matrix([[int(Aexp[x,p]-Aexp[0,p]) for p in range(r)] for x in range(1,M)])
    pr(f"{name}: M={M} r={r} rank of exponent differences = {diffs.rank()} (full = {M-1})")
    beta=(2*Aexp-1).T.astype(float)   # r x M, e_p=1
    n=M; d=n-1; T=1
    box=np.array(list(itertools.product(range(-T,T+1),repeat=d)))
    X=box@(beta[:,:d]-beta[:,[d]]).T
    mx=np.abs(X).max(axis=1); w=np.sort(mx)[int(0.5*len(mx))]
    A=box[mx<=w].astype(int); U=len(A)
    m=1
    while math.comb(m+1+n-2,n-1)*4<=U: m+=1
    mons=[a for t in range(m) for a in itertools.product(range(t+1),repeat=d) if sum(a)==t]
    Mx=sp.Matrix([[int(np.prod([int(k[j])**a[j] for j in range(d)])) for k in A] for a in mons])
    ns=Mx.nullspace()
    C=[]
    for v in ns[:60]:
        den=sp.ilcm(*[sp.fraction(x)[1] for x in v]); C.append([int(x*den) for x in v])
    H=(U*T**(m-1))**2
    pr(f"  slab-box: T={T} |A|={U} slab w={w} m={m} dim V(A,m)={len(ns)} (>= U/2: {len(ns)>=U/2}); "
       f"max sup-norm of first 60 basis vectors={max(max(abs(x) for x in c) for c in C)} <= H={H}")
    Afull=np.hstack([A,-A.sum(axis=1,keepdims=True)])
    viol=0; minrel=mp.inf
    Delta=mp.mpf('1e-6')
    for t in range(200):
        dl=[Delta*mp.mpf(float(rng.random())) for _ in range(M)]
        best=0
        for c in C:
            lam=[mp.fsum(int(kx)*dx for kx,dx in zip(k,dl)) for k in Afull]
            g=mp.fsum(ck*mp.e**(1j*l) for ck,l in zip(c,lam)) if False else None
            # leading Taylor term certifies nonvanishing: g = i^m Phi/m! + R, |R| <= sum|c| |lam|^{m+1}/(m+1)!
            Phi=mp.fsum(ck*l**m for ck,l in zip(c,lam))/mp.factorial(m)
            R=mp.fsum(abs(ck)*abs(l)**(m+1) for ck,l in zip(c,lam))/mp.factorial(m+1)
            if abs(Phi)>R:
                best=abs(Phi)/R; break
        if best>0: viol+=1; minrel=min(minrel,best)
    pr(f"  random fake positions delta~U[0,1e-6]^M: (4.2) certified violated (|leading term| > Taylor remainder) "
       f"for {viol}/200; min ratio leading/remainder = {mp.nstr(minrel,5)}")
open(sys.argv[1],"w").write("\n".join(out)+"\n")
