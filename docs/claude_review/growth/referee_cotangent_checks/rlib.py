# Independent exact library for the referee checks (Gaussian integers as int pairs).
from fractions import Fraction
from math import gcd, isqrt
import itertools, random

def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0],-a[1])
def gnorm(a): return a[0]*a[0]+a[1]*a[1]
def gdivmod_exact(a,b):
    # returns a/b if exact in Z[i], else None
    n=gnorm(b); num=gmul(a,gconj(b))
    if num[0]%n or num[1]%n: return None
    return (num[0]//n, num[1]//n)
def gmod(a,b):
    n=gnorm(b); num=gmul(a,gconj(b))
    q=( (2*num[0]+n)//(2*n), (2*num[1]+n)//(2*n) )
    r=(a[0]-gmul(q,b)[0], a[1]-gmul(q,b)[1]); return r
def ggcd(a,b):
    while b!=(0,0):
        a,b=b,gmod(a,b)
    return a
def ggcd_list(L):
    g=(0,0)
    for z in L: g=ggcd(g,z)
    return g

def factor(n):
    f={}; d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1
    if n>1: f[n]=f.get(n,0)+1
    return f

def gauss_prime_above(p):
    # p = 1 mod 4: find a+bi with a^2+b^2=p
    for a in range(1,isqrt(p)+1):
        b2=p-a*a; b=isqrt(b2)
        if b*b==b2: return (a,b)
    raise ValueError

def vpi(z,pi):
    v=0
    while True:
        q=gdivmod_exact(z,pi)
        if q is None: return v
        z=q; v+=1

def points_on_circle(N):
    pts=[]
    for x in range(-isqrt(N),isqrt(N)+1):
        y2=N-x*x; y=isqrt(y2)
        if y*y==y2:
            pts.append((x,y))
            if y: pts.append((x,-y))
    return pts

import math
def angle(z): return math.atan2(z[1],z[0])

def cot_half(zi,zj,N):
    # cot(theta_ij/2) = (N + <zi,zj>)/det(zi,zj)
    dot=zi[0]*zj[0]+zi[1]*zj[1]; det=zi[0]*zj[1]-zi[1]*zj[0]
    return Fraction(N+dot, det)

def primitive_cluster(pts):
    g=ggcd_list(pts)
    out=[gdivmod_exact(z,g) for z in pts]
    return out, g
