r"""Item (3)/(4d): heights of IX.3-type rational prime witnesses and the e-all singleton test at
POLYNOMIAL cutoffs.  Own code; nothing imported from Astra.
Construction (as in IX.3, with all marks of the class K_{i0}=1): Paley K_{ij}=chi(j-i) on F_q,
q=3 mod 4, labels I=F_q^*, a_{ij}=1+K_{ij}, marks S={i in I: K_{i0}=1}, anchor i0=min S.
D = rows a_i-a_{i0} (i in S\{i0}), H = rows e_i-e_{i0}.  Valid witness exponent matrices: integral V with
D V = (H-D)/2; F = I+2V.  Witness for label j: alpha_j = prod_l pi_l^{F_jl} (conjugates for negative
exponents), J_j = prod p_l^{|F_jl|}.  We pick V by an integral right inverse and then reduce each
column of V against the integral kernel of D (LLL + Babai), to make the witnesses as low as we can.
Then, with actual split primes p_l near q^d, we compute the odd X-smooth parts m_s of
Im a, Re a, Re a - Im a, Re a + Im a (the unit-1, -1, i, -i collision moduli of the singleton v=e_j)
and count labels with nu_s m_s > sqrt(p_j) (e-all fails for them whatever the angle)."""
import sys, math, random
from sympy import Matrix, ZZ, isprime, primerange
from sympy.polys.matrices import DomainMatrix
def chi(a,q):
    a%=q
    if a==0: return 0
    return 1 if pow(a,(q-1)//2,q)==1 else -1
def lll_split(D):
    m=len(D); n=len(D[0]); Wt=10**6
    rows=[[int(i==j) for j in range(n)]+[Wt*D[k][i] for k in range(m)] for i in range(n)]
    Bm=DomainMatrix([[ZZ(x) for x in r] for r in rows],(n,n+m),ZZ).lll().to_Matrix()
    ker=[];sol=[]
    for r in range(n):
        v=[int(x) for x in Bm.row(r)]
        (ker if all(x==0 for x in v[n:]) else sol).append(v)
    assert len(ker)==n-m and len(sol)==m
    Y=Matrix([[v[n+k]//Wt for v in sol] for k in range(m)])  # m x m, columns = D u_i
    assert abs(Y.det())==1
    Uc=Matrix([[v[i] for v in sol] for i in range(n)])        # n x m
    U=Uc*Y.inv()
    return U,[v[:n] for v in ker]
def right_inverse(D):
    return lll_split(D)[0]
def integral_kernel(D):
    return lll_split(D)[1]
def right_inverse_old(D):
    # column-style Hermite reduction with a unimodular transform T: D T = [L | 0]
    m=len(D); n=len(D[0])
    A=[row[:] for row in D]; T=[[int(i==j) for j in range(n)] for i in range(n)]
    def colop(c1,c2,a,b,c,d):  # (col c1, col c2) <- (a*c1+b*c2, c*c1+d*c2)
        for M_ in (A,T):
            for r in range(len(M_)):
                x,y=M_[r][c1],M_[r][c2]; M_[r][c1]=a*x+b*y; M_[r][c2]=c*x+d*y
    for r in range(m):
        for c in range(r+1,n):
            x,y=A[r][r],A[r][c]
            if y==0: continue
            # extended gcd
            g,s,t=egcd(x,y)
            colop(r,c,s,t,-y//g,x//g)
        assert abs(A[r][r])==1, ("invariant not one", r, A[r][r])
        if A[r][r]<0:
            for M_ in (A,T):
                for rr in range(len(M_)): M_[rr][r]*=-1
    # L = A[:, :m] lower triangular unit diagonal; U = T[:, :m] L^{-1}
    L=Matrix([row[:m] for row in A]); Linv=L.inv()
    U=Matrix([row[:m] for row in T])*Linv
    assert all(x==int(x) for x in U)
    return U
def egcd(a,b):
    if b==0: return (abs(a), 1 if a>=0 else -1, 0)
    g,s,t=egcd(b,a%b); return (g,t,s-(a//b)*t)
def gauss_prime(p):
    for x in range(1,math.isqrt(p)+1):
        y2=p-x*x; y=math.isqrt(y2)
        if y*y==y2: return (x,y)
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gpow(a,e):
    r=(1,0)
    while e:
        if e&1: r=gmul(r,a)
        a=gmul(a,a); e>>=1
    return r
def smooth_odd(n,X,excl):
    n=abs(n); s=1
    if n==0: return float('inf')
    for p in primerange(3,X+1):
        if p in excl: continue
        while n%p==0: n//=p; s*=p
    return s
def run(q,d=1.9,seed=0):
    random.seed(seed)
    I=list(range(1,q)); S=[i for i in I if chi(i-0,q)==1 and chi(0-i,q)==1] or [i for i in I if chi(0-i,q)==1]
    S=[i for i in I if chi(0-i,q)==1]  # K_{i0}=chi(0-i)
    i0=S[0]; m=len(S)
    a=lambda i,j: 1+chi(j-i,q)
    col={j:k for k,j in enumerate(I)}
    D=[[a(i,j)-a(i0,j) for j in I] for i in S[1:]]
    H=[[(1 if j==i else 0)-(1 if j==i0 else 0) for j in I] for i in S[1:]]
    import time; t0=time.time()
    U=right_inverse(D)
    print('U done',time.time()-t0,'max|U|',max(abs(int(x)) for x in U),flush=True)
    Dm=Matrix(D); Hm=Matrix(H)
    assert Dm*U==Matrix.eye(m-1)
    V=U*(Hm-Dm)/2
    assert all(x==int(x) for x in V)
    # kernel lattice of D (integral), LLL
    ker=Dm.nullspace()
    import functools
    kb=[]
    for v in ker:
        den=functools.reduce(lambda x,y: x*y//math.gcd(x,y), [int(x.q) for x in v],1)
        kb.append([int(x*den) for x in v])
    # saturate: integral kernel = {x in Z^n: Dx=0}; use HNF-free approach via T columns from right_inverse is cleaner
    # (columns m-1.. of T span the integral kernel); recompute:
    n=q-1
    K=integral_kernel(D)
    Kd=DomainMatrix([[ZZ(x) for x in row] for row in K],(len(K),n),ZZ).lll().to_Matrix()
    Kr=[list(map(int,Kd.row(r))) for r in range(Kd.rows)]
    F=Matrix.eye(n)+2*V
    # Babai rounding column by column on F columns: F col l + 2*k, k in kernel
    import numpy as np
    B=np.array(Kr,dtype=float)  # rows basis
    G=B@B.T; Ginv=np.linalg.inv(G)
    Fl=[[int(F[r,c]) for c in range(n)] for r in range(n)]
    assert Dm*Matrix(Fl)==Hm, 'before'
    assert all(all(x==0 for x in Dm*Matrix(k)) for k in Kr), 'kernel'
    for c in range(n):
        f=np.array([Fl[r][c] for r in range(n)],dtype=float)
        coef=np.rint(-(Ginv@(B@f))/2)
        ci=[int(x) for x in coef]; newf=[Fl[r][c]+2*sum(ci[k]*Kr[k][r] for k in range(len(Kr))) for r in range(n)]
        if sum(abs(x) for x in newf)<sum(abs(Fl[r][c]) for r in range(n)):
            for r in range(n): Fl[r][c]=newf[r]
    Fm=Matrix(Fl)
    assert Dm*Fm==Hm and all((Fm[r,c]-(r==c))%2==0 for r in range(n) for c in range(n))
    rown=[sum(abs(x) for x in Fl[r]) for r in range(n)]
    print('F done',time.time()-t0,'max|F|',max(abs(x) for r in Fl for x in r),'rown',sorted(rown),flush=True)
    # physical primes near P=q^d, distinct, split
    P=int(q**d); ps=[]; x=P
    while len(ps)<n:
        if isprime(x) and x%4==1: ps.append(x)
        x+=1
    random.shuffle(ps)
    pis=[gauss_prime(p) for p in ps]
    fails={}
    for X in (3*m, m*m):
        bad=0
        for j in range(n):
            al=(1,0)
            for l in range(n):
                e=Fl[j][l]
                if e>0: al=gmul(al,gpow(pis[l],e))
                elif e<0: al=gmul(al,gpow((pis[l][0],-pis[l][1]),-e))
            x_,y_=al
            ms=[smooth_odd(y_,X,set(ps)),smooth_odd(x_,X,set(ps)),smooth_odd(x_-y_,X,set(ps)),smooth_odd(x_+y_,X,set(ps))]
            nu=[1,1,2**-0.5,2**-0.5]
            if any(nu[s]*ms[s]>math.sqrt(ps[j]) for s in range(4)): bad+=1
        fails[X]=bad
    return m, n, max(abs(int(x)) for x in U), sorted(rown), fails, P
def integral_kernel(D):
    # integral kernel via the unimodular transform of the column reduction
    m=len(D); n=len(D[0])
    A=[row[:] for row in D]; T=[[int(i==j) for j in range(n)] for i in range(n)]
    for r in range(m):
        for c in range(r+1,n):
            x,y=A[r][r],A[r][c]
            if y==0: continue
            g,s,t=egcd(x,y)
            for M_ in (A,T):
                for rr in range(len(M_)):
                    u,v=M_[rr][r],M_[rr][c]; M_[rr][r]=s*u+t*v; M_[rr][c]=(-y//g)*u+(x//g)*v
    return [[T[r][c] for r in range(n)] for c in range(m,n)]
for q in [int(a) for a in sys.argv[1:]] or [11,19,23,31,43]:
    m,n,umax,rown,fails,P=run(q)
    print(f"q={q} m={m} labels={n} P~{P}  max|U|={umax}  row L1 norms of F (log J_j/log p): min {rown[0]} median {rown[len(rown)//2]} max {rown[-1]}  "
          f"labels failing singleton e-all for every angle: X=3m:{fails[3*m]}/{n}, X=m^2:{fails[m*m]}/{n}")
