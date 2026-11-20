# Referee check of the X.5 negative control (labels 101, 109, prime 1511).
# Exact arithmetic in F_l[i] = Z[i]/l, l = 3 mod 4 (inert) or split l.
from sympy import isprime, primerange, sqrt_mod
from math import isqrt
import itertools, random

def mul(a, b, l):
    return ((a[0]*b[0]-a[1]*b[1]) % l, (a[0]*b[1]+a[1]*b[0]) % l)
def conj(a, l): return (a[0] % l, (-a[1]) % l)
def nrm(a, l): return (a[0]*a[0]+a[1]*a[1]) % l
def inv(a, l):
    n = nrm(a, l); ni = pow(n, -1, l)
    c = conj(a, l); return (c[0]*ni % l, c[1]*ni % l)
def ratio(a, l): return mul(a, inv(conj(a, l), l), l)
def leg(a, l): return pow(a % l, (l-1)//2, l)

# 1. the prime 1511
l = 1511
assert isprime(l) and l % 4 == 3
primes = [p for p in primerange(101, 2000) if p % 4 == 1][:80]
print("80 fine labels from", primes[0], "to", primes[-1], "; 1511 parent?", l in primes)
first = next(e for e in primerange(1500, 4500) if leg(101, e) == leg(109, e))
print("first prime >=1500 with equal Legendre symbols of 101,109:", first,
      "symbols", leg(101, l), leg(109, l))
# 2. exhaustive: all units of norm 101 and of norm 109 mod 1511; ratio sets
R1 = {}; R2 = {}
for x in range(l):
    for y in range(l):
        n = (x*x+y*y) % l
        if n == 101 % l: R1.setdefault(ratio((x, y), l), (x, y))
        elif n == 109 % l: R2.setdefault(ratio((x, y), l), (x, y))
common = set(R1) & set(R2)
print("ratio classes norm101:", len(R1), " norm109:", len(R2), " common:", len(common),
      " (norm-one group order l+1 =", l+1, ")")
r1, r2 = R1[min(common)], R2[min(common)]
q = mul(ratio(r1, l), inv(ratio(r2, l), l), l)
print("explicit r_101 =", r1, " r_109 =", r2, " (r/rbar)_101 / (r/rbar)_109 =", q)
assert q == (1, 0)
gm = 101*109
print("sqrt(101*109) =", gm**0.5, "; 1511^2 > 101*109:", l*l > gm,
      "; required |sin| >= 1511/sqrt(11009) =", l/gm**0.5)

# 3. criterion: equal ratios possible iff equal Legendre symbols, at small inert and split primes
def ratioset(p, l):
    return {ratio((x, y), l) for x in range(l) for y in range(l) if (x*x+y*y) % l == p % l}
bad = 0
for l2 in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]:
    for p, pp in itertools.combinations([5, 13, 17, 29, 37, 41, 53, 61], 2):
        if l2 in (p, pp): continue
        inter = bool(ratioset(p, l2) & ratioset(pp, l2))
        if inter != (leg(p, l2) == leg(pp, l2)): bad += 1
print("Legendre criterion for equal ratios, small-prime exhaustive mismatches:", bad)
