#!/usr/bin/env python3
"""Pairs z1=g(a+i), z2=g(a-i), |g|^2 = n prime = 1 mod 4 just below a^2: angular distance vs 2 N^{-1/4},
primitive height H of the pair, and log H / log prod_{p|N} p (N squarefree)."""
import math
from sympy import prevprime, factorint
def gp(p):
    for y in range(1, math.isqrt(p)+1):
        x = math.isqrt(p-y*y)
        if x*x+y*y == p: return (x, y)
def mul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
for a in (10**3, 10**5, 10**7, 10**9):
    n = prevprime(a*a)
    while n % 4 != 1: n = prevprime(n)
    g = gp(n); N = n*(a*a+1)
    z1 = mul(g, (a, 1)); z2 = mul(g, (a, -1))
    ang = abs(math.atan2(z1[1], z1[0]) - math.atan2(z2[1], z2[0]))
    fac = factorint(a*a+1); sqf = all(e == 1 for e in fac.values())
    rad = n*math.prod(fac)
    H = a*a+1   # gcd(z1,z2)=g since gcd(a+i,a-i)=1 for a even
    print('a=%g  C_needed=ang*N^(1/4)=%.4f  a^2+1 squarefree=%s  logH/log rad(N)=%.6f' % (a, ang*N**0.25, sqf, math.log(H)/math.log(rad)))
