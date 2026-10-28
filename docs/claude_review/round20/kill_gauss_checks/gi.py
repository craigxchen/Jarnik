# minimal exact Gaussian-integer toolkit (own code)
from math import isqrt
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0],-a[1])
def gnorm(a): return a[0]*a[0]+a[1]*a[1]
def rdiv(n,d):  # round(n/d) for ints
    q,r=divmod(2*n+d,2*d); return q
def gdivmod(a,b):
    nb=gnorm(b); num=gmul(a,gconj(b))
    q=(rdiv(num[0],nb), rdiv(num[1],nb))
    qb=gmul(q,b); return q,(a[0]-qb[0],a[1]-qb[1])
def gexact(a,b):
    q,r=gdivmod(a,b); assert r==(0,0),(a,b,r); return q
def ggcd(a,b):
    while b!=(0,0):
        _,r=gdivmod(a,b); a,b=b,r
    return a
def ggcdlist(L):
    g=(0,0)
    for z in L: g=ggcd(g,z) if g!=(0,0) else z
    return g
