# Exact checks for Lemma 5.2 and Lemma 5.3 (analytic.md section 5).
# (1) Lemma 5.2: for every Gaussian prime pi=a+bi (a>b>0) with p = a^2+b^2 <= 2000 and every alpha in Z[i],
#     0 < Norm alpha <= 3000, not parallel to pi: |Im(pi conj alpha)| >= 1 (integer), i.e. LCG_p holds.
# (2) Lemma 5.3 (contrapositive, exact rational arithmetic): for directions phi = arg(u+iv) with u+iv primitive and
#     Norm(u+iv) > n, there is alpha with Norm alpha <= n+1, not parallel, violating LCG_n; for Norm(u+iv) <= n none.
#     (sin(phi-arg alpha)^2 = Im(alpha conj w)^2/(Norm alpha Norm w), compared exactly with 1/(n Norm alpha).)
from fractions import Fraction
from math import gcd, isqrt
from sympy import isprime
viol = 0; tested = 0
alphas = [(x, y) for x in range(-55, 56) for y in range(-55, 56) if 0 < x*x + y*y <= 3000]
for a in range(1, 45):
    for b in range(1, a):
        p = a*a + b*b
        if p > 2000 or not isprime(p): continue
        for (x, y) in alphas:
            im = b*x - a*y          # Im((a+bi)(x-iy)) = b x - a y
            if im == 0: continue   # parallel
            tested += 1
            if im*im < 1: viol += 1
print("Lemma 5.2: pairs tested %d, violations %d" % (tested, viol))

def lcg_violator(w, n):
    u, v = w; Nw = u*u + v*v
    best = None
    R = isqrt(n + 1)
    for x in range(-R - 1, R + 2):
        for y in range(-R - 1, R + 2):
            Na = x*x + y*y
            if Na == 0 or Na > n + 1: continue
            im = v*x - u*y         # Im(w conj alpha) up to sign
            if im == 0: continue
            # violation iff im^2/(Na*Nw) < 1/(n*Na)  iff  n*im^2 < Nw
            if n*im*im < Nw: return (x, y)
    return None
bad = 0; cnt = 0
for n in (5, 13, 17, 29, 101):
    for u in range(-40, 41):
        for v in range(-40, 41):
            if (u, v) == (0, 0) or gcd(abs(u), abs(v)) != 1: continue
            cnt += 1
            viol = lcg_violator((u, v), n)
            Nw = u*u + v*v
            if Nw > n and viol is None: bad += 1; print("no violator for", (u, v), n)
            if Nw <= n and viol is not None: bad += 1; print("unexpected violator", (u, v), n, viol)
print("Lemma 5.3: %d (direction, n) pairs, exceptions %d" % (cnt, bad))
