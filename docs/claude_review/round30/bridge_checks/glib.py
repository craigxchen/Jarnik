# exact helpers: Gaussian rationals (pairs of Fractions), polynomials over Q(i) as lists (low->high), Z[i] gcd.
from fractions import Fraction as Fr
from math import gcd
import itertools

class G:
    __slots__=('a','b')
    def __init__(s,a,b=0): s.a=Fr(a); s.b=Fr(b)
    def __add__(s,o): o=g_(o); return G(s.a+o.a,s.b+o.b)
    __radd__=__add__
    def __sub__(s,o): o=g_(o); return G(s.a-o.a,s.b-o.b)
    def __rsub__(s,o): return g_(o)-s
    def __neg__(s): return G(-s.a,-s.b)
    def __mul__(s,o): o=g_(o); return G(s.a*o.a-s.b*o.b,s.a*o.b+s.b*o.a)
    __rmul__=__mul__
    def conj(s): return G(s.a,-s.b)
    def norm(s): return s.a*s.a+s.b*s.b
    def __truediv__(s,o):
        o=g_(o); n=o.norm(); c=s*o.conj(); return G(c.a/n,c.b/n)
    def __eq__(s,o): o=g_(o); return s.a==o.a and s.b==o.b
    def __hash__(s): return hash((s.a,s.b))
    def iszero(s): return s.a==0 and s.b==0
    def __repr__(s): return f"({s.a}+{s.b}i)"
def g_(x): return x if isinstance(x,G) else G(x,0)

def pmul(p,q):
    r=[G(0)]*(len(p)+len(q)-1)
    for i,x in enumerate(p):
        for j,y in enumerate(q): r[i+j]=r[i+j]+x*y
    return r
def psub(p,q):
    n=max(len(p),len(q)); p=p+[G(0)]*(n-len(p)); q=q+[G(0)]*(n-len(q))
    return trim([x-y for x,y in zip(p,q)])
def trim(p):
    p=list(p)
    while len(p)>1 and p[-1].iszero(): p.pop()
    return p
def deg(p):
    p=trim(p); return -1 if (len(p)==1 and p[0].iszero()) else len(p)-1
def pconj(p): return [x.conj() for x in p]
def peval(p,t):
    t=g_(t); r=G(0)
    for c in reversed(p): r=r*t+c
    return r
def prod_lin(roots):
    p=[G(1)]
    for r in roots: p=pmul(p,[-g_(r),G(1)])
    return p

# Z[i] arithmetic on python int pairs
def gi_divmod(a,b):
    # a,b tuples (x,y); nearest-integer quotient
    n=b[0]*b[0]+b[1]*b[1]
    num=(a[0]*b[0]+a[1]*b[1], a[1]*b[0]-a[0]*b[1])
    def rnd(x): return (2*x+n)//(2*n)
    q=(rnd(num[0]),rnd(num[1]))
    r=(a[0]-(q[0]*b[0]-q[1]*b[1]), a[1]-(q[0]*b[1]+q[1]*b[0]))
    return q,r
def gi_gcd(a,b):
    while b!=(0,0):
        _,r=gi_divmod(a,b); a,b=b,r
    return a
def frac_gcd(vals):
    """vals: list of G (Gaussian rationals). Return (num tuple, den int) with ideal gcd = num/den."""
    den=1
    for v in vals:
        for q in (v.a,v.b): den=den*q.denominator//gcd(den,q.denominator)
    g=(0,0)
    for v in vals:
        x=(int(v.a*den),int(v.b*den)); g=gi_gcd(g,x) if g!=(0,0) else x
    return g,den
def gnorm_int(t): return t[0]*t[0]+t[1]*t[1]
