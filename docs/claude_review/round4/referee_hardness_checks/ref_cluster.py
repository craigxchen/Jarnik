#!/usr/bin/env python3
"""Referee check R1: the 10-point cluster on N = 1176852625, exact Gaussian arithmetic.
 (a) list all 45 pairs sorted by ANGLE (not by pair norm), with pair norm n, |gcd|, ratio to sqrt2/sqrt n;
 (b) residue classes: for each pair, the gcd g gives a modulus with z = w (mod g).  Report the
     exponent alpha = log Norm(g) / log Norm(N) = log|g|^2/(2 log N) and compare with 1/4 and with
     1/4 + log C/(2 log N)... (Prop 6.1 needs |g| > C N^(1/4)).
 (c) exact: arc constant of the 10 points.
"""
import math, itertools
from fractions import Fraction

def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0], -a[1])
def gnorm(a): return a[0]*a[0]+a[1]*a[1]
def gdivexact(a, b):
    nb = gnorm(b); num = gmul(a, gconj(b))
    assert num[0] % nb == 0 and num[1] % nb == 0
    return (num[0]//nb, num[1]//nb)
def gdivmod_round(a, b):
    nb = gnorm(b); num = gmul(a, gconj(b))
    q = ((2*num[0] + nb)//(2*nb), (2*num[1] + nb)//(2*nb))
    qb = gmul(q, b)
    return q, (a[0]-qb[0], a[1]-qb[1])
def ggcd(a, b):
    while b != (0, 0):
        _, r = gdivmod_round(a, b); a, b = b, r
    return a

N = 1176852625
r = math.isqrt(N)
pts = []
for x in range(-r, r+1):
    y2 = N - x*x; y = math.isqrt(y2)
    if y*y == y2:
        pts.append((x, y))
        if y: pts.append((x, -y))
R = math.sqrt(N)
ang = sorted((math.atan2(p[1], p[0]) % (2*math.pi), p) for p in pts)
m = len(ang)
best = None
for i in range(m):
    j = i + 9
    a1 = ang[j % m][0] + (2*math.pi if j >= m else 0)
    arc = R*(a1 - ang[i][0])
    if best is None or arc < best[0]:
        best = (arc, [ang[(i+t) % m][1] for t in range(10)])
C = best[0]/math.sqrt(R)
P = best[1]
print(f"N={N}, r2(N)={len(pts)}, 10-point arc constant C={C:.6f}")
rows = []
for z, w in itertools.combinations(P, 2):
    g = ggcd(z, w)
    n = N // gnorm(g); assert N % gnorm(g) == 0
    q = gmul(z, gconj(w))
    lam = abs(math.atan2(q[1], q[0]))
    # exact check that z = w mod g
    d = (z[0]-w[0], z[1]-w[1]); gdivexact(d, g)
    alpha = math.log(gnorm(g))/(2*math.log(N))
    rows.append((lam, n, gnorm(g), lam/(math.sqrt(2)/math.sqrt(n)), alpha, abs(complex(*d))))
rows.sort()
print("eight smallest angles (= eight closest pairs):")
print("  angle      pair-norm n   Norm(g)     ratio   alpha=logNorm(g)/logNorm(N)  |g|/(C N^(1/4))")
for lam, n, ng, ratio, alpha, dist in rows[:8]:
    print(f"  {lam:.5e} {n:>12} {ng:>10} {ratio:8.3f}   {alpha:.4f}   {math.sqrt(ng)/(C*N**0.25):.3f}")
print("pairs sorted by pair norm, first 8 (what check_forms.py printed):")
for lam, n, ng, ratio, alpha, dist in sorted(rows, key=lambda t: (t[1], t[0]))[:8]:
    print(f"  n={n:>10} angle={lam:.5e} ratio={ratio:.3f}")
cnt = sum(1 for t in rows if t[4] > 0.25)
print(f"pairs whose gcd modulus g has alpha > 1/4 (so two cluster points share the class 0 mod g "
      f"with Norm(g) > Norm(N)^(1/4)): {cnt} of 45")
mx = max(rows, key=lambda t: t[4])
print(f"  largest such alpha = {mx[4]:.4f} (Norm g = {mx[2]}), versus 1/4 + log C/log N = "
      f"{0.25 + math.log(C)/math.log(N):.4f}")
