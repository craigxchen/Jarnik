"""Judge check J1 (exact). Sylvester profile on F_2^k, columns beta_u(x) = (-1)^{u.x}, u != 0, uniform weights.
(a) every x -> Ax+b (A in GL(k,2)) maps each column to +-another column (all of AGL(3,2); generators + 3000 random of AGL(4,2));
(b) Lambda is Aut-invariant: Lambda_sigma(c) = Lambda(c) for random supports S (h_S(v) = max_S<k,v> - min_S<k,v>);
(c) identity P_c(Z^sigma) = P_{c o sigma}(Z) in exact Gaussian rationals (the certificate at Z^sigma is the
    route-5 certificate at Z for the permuted form);
(d) exact ln|AGL(12,2)|, ln|AGL(13,2)|."""
import itertools, random, math
from fractions import Fraction as Fr
random.seed(26)
def pts(k): return list(itertools.product((0,1),repeat=k))
def dot(u,x): return sum(a*b for a,b in zip(u,x))%2
def col(u,P): return tuple((-1)**dot(u,x) for x in P)
def apply(A,b,x): return tuple((sum(A[i][j]*x[j] for j in range(len(x)))+b[i])%2 for i in range(len(x)))
def rank2(A):
    M=[r[:] for r in A]; r=0; k=len(M)
    for c in range(k):
        piv=next((i for i in range(r,k) if M[i][c]),None)
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]
        for i in range(k):
            if i!=r and M[i][c]: M[i]=[(a+b)%2 for a,b in zip(M[i],M[r])]
        r+=1
    return r
def check_sigma(k,A,b,P,cols):
    idx={x:i for i,x in enumerate(P)}
    perm=[idx[apply(A,b,x)] for x in P]          # sigma: point i -> perm[i]
    colset=set(cols)|set(tuple(-v for v in c) for c in cols)
    for c in cols:
        img=[0]*len(P)
        for i in range(len(P)): img[perm[i]]=c[i]
        if tuple(img) not in colset: return False,perm
    return True,perm
def h(S,v): vals=[sum(a*b for a,b in zip(kk,v)) for kk in S]; return max(vals)-min(vals)
for k in (3,4):
    P=pts(k); n=len(P); U=[u for u in P if any(u)]; cols=[col(u,P) for u in U]
    mats=[]
    if k==3:
        for flat in itertools.product((0,1),repeat=9):
            A=[list(flat[3*i:3*i+3]) for i in range(3)]
            if rank2(A)==3: mats.append(A)
        group=[(A,b) for A in mats for b in P]
    else:
        group=[]
        while len(group)<3000:
            A=[[random.randint(0,1) for _ in range(k)] for _ in range(k)]
            if rank2(A)==k: group.append((A,random.choice(P)))
    ok=0; perms=[]
    for A,b in group:
        good,perm=check_sigma(k,A,b,P,cols); ok+=good; perms.append(perm)
    print(f"k={k} n={n}: checked {len(group)} affine maps, column-to-+-column for {ok}")
    # (b) Lambda invariance
    bad=0
    for t in range(200):
        S=[]
        while len(S)<5:
            kk=[random.randint(-2,2) for _ in range(n-1)]; kk.append(-sum(kk)); S.append(tuple(kk))
        perm=random.choice(perms)
        Ssig=[tuple(kk[perm.index(i)] for i in range(n)) for kk in S]   # k o sigma^{-1}
        L1=sum(h(S,c) for c in cols); L2=sum(h(Ssig,c) for c in cols)
        bad+= (L1!=L2)
    print(f"  Lambda(sigma S) == Lambda(S) on 200 random supports: failures={bad}")
    # (c) identity P_c(Z^sigma) = P_{c o sigma}(Z), exact Gaussian rationals as pairs of Fractions
    def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
    def ginv(a): N=a[0]**2+a[1]**2; return (a[0]/N, -a[1]/N)
    def gpow(a,e):
        r=(Fr(1),Fr(0)); base=a if e>=0 else ginv(a)
        for _ in range(abs(e)): r=gmul(r,base)
        return r
    def mono(Z,kk):
        r=(Fr(1),Fr(0))
        for z,e in zip(Z,kk): r=gmul(r,gpow(z,e))
        return r
    bad=0
    for t in range(30):
        Z=[(Fr(random.randint(1,9)),Fr(random.randint(-9,9))) for _ in range(n)]
        perm=random.choice(perms)
        Zs=[Z[perm[i]] for i in range(n)]                  # (Z^sigma)_i = z_{sigma(i)}
        C={}
        for _ in range(4):
            kk=[random.randint(-2,2) for _ in range(n-1)]; kk.append(-sum(kk)); C[tuple(kk)]=random.randint(-5,5)
        lhs=(Fr(0),Fr(0))
        for kk,cv in C.items():
            m=mono(Zs,kk); lhs=(lhs[0]+cv*m[0],lhs[1]+cv*m[1])
        rhs=(Fr(0),Fr(0))
        for kk,cv in C.items():
            kp=[0]*n
            for i in range(n): kp[perm[i]]=kk[i]           # exponent moved to coordinate sigma(i)
            m=mono(Z,kp); rhs=(rhs[0]+cv*m[0],rhs[1]+cv*m[1])
        bad+=(lhs!=rhs)
    print(f"  P_c(Z^sigma) == P_(c o sigma^-1 on exponents)(Z), 30 exact trials: failures={bad}")
def agl(k):
    g=2**k
    for i in range(k): g*= (2**k-2**i)
    return g
for k in (12,13): print(f"ln|AGL({k},2)| = {math.log(agl(k)):.2f}")
