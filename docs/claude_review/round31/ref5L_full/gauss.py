"""Small Gaussian-integer toolkit (independent of Astra's scripts)."""
from math import gcd
def mul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def conj(a): return (a[0],-a[1])
def norm(a): return a[0]*a[0]+a[1]*a[1]
def rdiv(n,d):  # nearest integer to n/d, d>0
    return (2*n+d)//(2*d)
def divmod_g(a,b):
    nb=norm(b); num=mul(a,conj(b))
    q=(rdiv(num[0],nb), rdiv(num[1],nb))
    r=(a[0]-mul(q,b)[0], a[1]-mul(q,b)[1])
    return q,r
def ggcd(a,b):
    while b!=(0,0):
        _,r=divmod_g(a,b); a,b=b,r
    return a
def exdiv(a,b):
    q,r=divmod_g(a,b); assert r==(0,0); return q
def divides(b,a):
    return divmod_g(a,b)[1]==(0,0)
def det(a,b): return a[0]*b[1]-a[1]*b[0]
def dot(a,b): return a[0]*b[0]+a[1]*b[1]
def genuine(Vs):
    """Least-norm primitive Gaussian points zeta_i = lam V_i/conj(V_i). Returns (zetas,H)."""
    Lam=1
    for v in Vs: Lam*=norm(v)
    Z=[]
    for v in Vs:
        s=mul(v,v); f=Lam//norm(v); Z.append((s[0]*f,s[1]*f))
    g=Z[0]
    for z in Z[1:]: g=ggcd(g,z)
    zs=[exdiv(z,g) for z in Z]
    H=norm(zs[0]); assert all(norm(z)==H for z in zs)
    return zs,H
