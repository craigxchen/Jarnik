#!/usr/bin/env python3
"""Independent brute force: enumerate b = x+iy directly by x,y with h=x^2+y^2 <= B, P-supported,
gcd(b, conj b)=1, and test the integer form directly.  Counts pairs (b, eps) / 4 = number of u."""
import sys, math
from fractions import Fraction
from sympy import primerange
s = int(sys.argv[1]); B = int(sys.argv[2]); delta = Fraction(sys.argv[3]); gamma = Fraction(sys.argv[4])
ps = [p for p in primerange(5, 10**6) if p % 4 == 1][:s]
prodP = math.prod(ps)
def supported(h):
    for p in ps:
        while h % p == 0: h //= p
    return h == 1
cnt = 0
R = math.isqrt(B)
for x in range(-R, R+1):
    for y in range(-R, R+1):
        h = x*x+y*y
        if h == 0 or h > B: continue
        if h**gamma.denominator < prodP**gamma.numerator: continue
        if not supported(h): continue
        if math.gcd(x, y) != 1 or (x - y) % 2 == 0: continue   # gcd(b, conj b)=1 <=> gcd(x,y)=1 and x+y odd
        for L in (4*y*y, 4*x*x, 2*(x+y)**2, 2*(x-y)**2):
            if L > 0 and L**delta.denominator <= h**(delta.denominator-delta.numerator):
                cnt += 1
print('brute: pairs(b,eps)=%d  u-count=%s' % (cnt, Fraction(cnt, 4)))
