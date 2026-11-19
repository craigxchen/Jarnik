"""Exact Gaussian integer helpers (pairs of Python ints)."""
from math import isqrt
def mul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def conj(a): return (a[0],-a[1])
def nrm(a): return a[0]*a[0]+a[1]*a[1]
def rdiv(n,d):  # nearest integer to n/d, d>0
    return (2*n+d)//(2*d)
def divmod_g(a,b):
    n=nrm(b); t=mul(a,conj(b)); q=(rdiv(t[0],n),rdiv(t[1],n))
    p=mul(q,b); return q,(a[0]-p[0],a[1]-p[1])
def gcd(a,b):
    while b!=(0,0):
        _,r=divmod_g(a,b); a,b=b,r
    return a
def exdiv(a,b):
    q,r=divmod_g(a,b); assert r==(0,0),(a,b); return q
def lcm(a,b):
    g=gcd(a,b); return exdiv(mul(a,b),g)
def split_prime(p):
    """pi with N(pi)=p for p=2 or p=1 mod 4"""
    if p==2: return (1,1)
    # find sqrt(-1) mod p
    for c in range(2,p):
        if pow(c,(p-1)//2,p)==p-1:
            r=pow(c,(p-1)//4,p); break
    return gcd((p,0),(r,1))
def val(a,pi):
    v=0
    while True:
        q,r=divmod_g(a,pi)
        if r!=(0,0): return v
        a=q; v+=1
