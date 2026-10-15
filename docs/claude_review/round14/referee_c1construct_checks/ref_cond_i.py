# Referee check: how often does condition (i) of Prop. P fail on prime sets a fake may choose?
# Independent code (no import from the author's scripts).
from sympy import primerange, isprime
from math import isqrt, atan2, pi, log
import random

def gauss(p):
    # a odd >0, b even >0, a^2+b^2=p
    for b in range(2, isqrt(p)+1, 2):
        a2 = p - b*b
        a = isqrt(a2)
        if a*a == a2 and a % 2 == 1:
            return a, b
    raise ValueError(p)

def S_of(a, b):
    th = atan2(b, a)
    assert abs(th - pi/4) > 1e-12
    return (2, 3) if th < pi/4 else (1, 2, 3)

def D(s, a, b):
    return {0: b, 1: a-b, 2: a, 3: a+b}[s]

def fails_i(j, ps, gs):
    a, b = gs[j]
    prod = 1
    for s in S_of(a, b):
        prod *= abs(D(s, a, b))
    return any(prod % pl == 0 for pl in ps)

def count_fail(ps):
    gs = [gauss(p) for p in ps]
    return sum(fails_i(j, ps, gs) for j in range(len(ps)))

split = [p for p in primerange(3, 3*10**6) if p % 4 == 1]
print("(A) first r split primes (the near-extremal prime sets of actual circles):")
for r in [10, 20, 50, 100, 200, 500, 1000]:
    ps = split[:r]
    print(f"   r={r:5d}, p_max={ps[-1]:7d}: columns failing (i) = {count_fail(ps)} / {r}")

print("(B) selected sets: 5 and 13 plus r-2 primes p<=Y chosen with 5 or 13 dividing prod_{s in S_j} D_s(pi_j):")
for r, Y in [(50, 10**4), (200, 10**5), (1000, 10**6)]:
    base = [5, 13]
    cand = [p for p in split if p > 13 and p <= Y]
    good = []
    for p in cand:
        a, b = gauss(p)
        pr = 1
        for s in S_of(a, b):
            pr *= abs(D(s, a, b))
        if pr % 5 == 0 or pr % 13 == 0:
            good.append(p)
    frac = len(good)/len(cand)
    random.seed(1)
    ps = base + random.sample(good, r-2)
    print(f"   Y={Y}: fraction of split primes <= Y with 5|.. or 13|..: {frac:.3f}; chosen r={r}: failing (i) = {count_fail(ps)} / {r}")
