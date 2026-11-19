# Exact Gaussian-integer utilities, written independently for the referee check.
from math import gcd, isqrt
from fractions import Fraction
from itertools import combinations
import random

def gmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0],-a[1])
def gnorm(a): return a[0]*a[0]+a[1]*a[1]
def gpow(a,k):
    r=(1,0)
    for _ in range(k): r=gmul(r,a)
    return r
def gdivmod_exact(a,b):
    # return a/b if exact in Z[i], else None
    n=gnorm(b); num=gmul(a,gconj(b))
    if num[0]%n or num[1]%n: return None
    return (num[0]//n, num[1]//n)
def ground_div(a,b):
    n=gnorm(b); num=gmul(a,gconj(b))
    def rd(x): return (2*x+n)//(2*n)
    q=(rd(num[0]),rd(num[1]))
    r=(a[0]-gmul(q,b)[0], a[1]-gmul(q,b)[1])
    return q,r
def ggcd(a,b):
    while b!=(0,0):
        q,r=ground_div(a,b); a,b=b,r
    return a
def val(z,pi):
    k=0
    while True:
        q=gdivmod_exact(z,pi)
        if q is None: return k
        z=q; k+=1
def split_prime(p):
    # find pi=a+bi with a^2+b^2=p
    for a in range(1,isqrt(p)+1):
        b2=p-a*a; b=isqrt(b2)
        if b*b==b2: return (a,b)
SPLIT=[p for p in range(5,400) if p%4==1 and all(p%q for q in range(2,isqrt(p)+1))]
PI={p:split_prime(p) for p in SPLIT}
UNITS=[(1,0),(0,1),(-1,0),(0,-1)]
