# Referee check of Lemma P0 (independent): for theta = arg(a+ib), the arcs |sin(x - s pi/4)| < |sin(theta - s pi/4)|,
# s in S(theta), cover [0,pi) except finitely many points; list the uncovered points and compare with the claim.
from sympy import primerange
from math import isqrt
import mpmath as mp
mp.mp.dps = 50
def gauss(p):
    for b in range(2, isqrt(p)+1, 2):
        a2 = p-b*b; a = isqrt(a2)
        if a*a == a2 and a % 2: return a, b
bad = 0; n = 0
for p in primerange(5, 10**5):
    if p % 4 != 1: continue
    a, b = gauss(p); th = mp.atan2(b, a)
    S = (2, 3) if th < mp.pi/4 else (1, 2, 3)
    ivs = []
    for s in S:
        c = s*mp.pi/4; r = abs(mp.asin(mp.sin(th - c)))   # dist(theta - c, pi Z)
        for sh in (-mp.pi, 0, mp.pi):
            L, H = max(c-r+sh, mp.mpf(0)), min(c+r+sh, mp.pi)
            if L < H: ivs.append((L, H))
    ivs.sort(); cur = mp.mpf(0); free = []
    for L, H in ivs:
        if L > cur + mp.mpf(10)**-40: free.append((cur, L))
        cur = max(cur, H)
    if cur < mp.pi - mp.mpf(10)**-40: free.append((cur, mp.pi))
    # endpoints (zero-length free points) where arcs touch: theta and possibly pi/2-theta
    n += 1
    if free:   # positive-length free interval => lemma false
        bad += 1; print("POSITIVE-LENGTH FREE SET", p, free)
print(f"Lemma P0: {n} split primes < 1e5 checked; positive-length free sets: {bad}")
