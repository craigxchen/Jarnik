# Own exact arithmetic over Q(i)[t]: polynomials as lists of Gaussian rationals (pairs of Fractions)
from fractions import Fraction as Fr
import itertools, math

class G:  # Gaussian rational
    __slots__=('a','b')
    def __init__(s,a,b=0): s.a=Fr(a); s.b=Fr(b)
    def __add__(s,o): return G(s.a+o.a,s.b+o.b)
    def __sub__(s,o): return G(s.a-o.a,s.b-o.b)
    def __mul__(s,o): return G(s.a*o.a-s.b*o.b,s.a*o.b+s.b*o.a)
    def __neg__(s): return G(-s.a,-s.b)
    def conj(s): return G(s.a,-s.b)
    def iszero(s): return s.a==0 and s.b==0
    def inv(s):
        d=s.a*s.a+s.b*s.b; return G(s.a/d,-s.b/d)
    def __eq__(s,o): return s.a==o.a and s.b==o.b
    def __repr__(s): return f"({s.a}+{s.b}i)"
Z0=G(0); ONE=G(1); I=G(0,1)

def trim(p):
    p=list(p)
    while p and p[-1].iszero(): p.pop()
    return p
def padd(p,q):
    n=max(len(p),len(q)); return trim([(p[k] if k<len(p) else Z0)+(q[k] if k<len(q) else Z0) for k in range(n)])
def psub(p,q): return padd(p,[-c for c in q])
def pmul(p,q):
    if not p or not q: return []
    r=[Z0]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        if a.iszero(): continue
        for j,b in enumerate(q): r[i+j]=r[i+j]+a*b
    return trim(r)
def pconj(p): return [c.conj() for c in p]
def deg(p): return len(p)-1  # deg(0) = -1
def pdivmod(p,q):
    p=trim(p); q=trim(q); qi=q[-1].inv(); quo=[Z0]*max(0,len(p)-len(q)+1)
    p=list(p)
    while len(p)>=len(q) and p:
        c=p[-1]*qi; k=len(p)-len(q); quo[k]=c
        for j,b in enumerate(q): p[k+j]=p[k+j]-c*b
        p=trim(p)
    return trim(quo),p
def pgcd(p,q):
    p=trim(p); q=trim(q)
    while q:
        _,r=pdivmod(p,q); p,q=q,r
    if p: inv=p[-1].inv(); p=[c*inv for c in p]
    return p
def gcdlist(L):
    g=L[0]
    for x in L[1:]: g=pgcd(g,x)
    return g
def P(*coeffs): return trim([c if isinstance(c,G) else G(c) for c in coeffs])  # low to high
T=P(0,1)
def matmul(A,B):
    return [[padd(pmul(A[i][0],B[0][j]),pmul(A[i][1],B[1][j])) for j in range(2)] for i in range(2)]
def matpow(A,n):
    R=[[P(1),[]],[[],P(1)]]
    for _ in range(n): R=matmul(R,A)
    return R
def block(An,s):
    s=G(s)
    X=padd(An[0][0],pmul(An[0][1],[s])); Y=padd(An[1][0],pmul(An[1][1],[s]))
    return padd(X,pmul(Y,[I]))
def prod(L):
    r=P(1)
    for x in L: r=pmul(r,x)
    return r
def analyse(marks):
    g=gcdlist(marks)
    red=[pdivmod(z,g)[0] for z in marks]
    norms=[pmul(z,pconj(z)) for z in marks]
    eqn=all(norms[0]==nn for nn in norms)
    # normalise so that leading coeffs equal (projective tuple): contact defined via ratios
    D=deg(red[0]); assert all(deg(z)==D for z in red)
    lc=[z[-1] for z in red]
    nz=[pmul(z,[lc[0]*lc[k].inv()]) for k,z in enumerate(red)]
    md=-1; pairc={}
    for i,j in itertools.combinations(range(len(nz)),2):
        d=deg(psub(nz[i],nz[j])); pairc[(i,j)]=D-d; md=max(md,d)
    distinct=all(psub(nz[i],nz[j]) for i,j in itertools.combinations(range(len(nz)),2))
    return dict(eqnorm=eqn,gcddeg=deg(g),D=D,h=D-md,distinct=distinct,pairc=pairc,lcmod=[ (c.a*c.a+c.b*c.b) for c in lc])
