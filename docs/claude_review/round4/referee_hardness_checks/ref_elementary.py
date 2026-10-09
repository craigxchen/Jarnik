#!/usr/bin/env python3
"""Referee checks R2-R7 (exact integer arithmetic unless stated).
R2 Lemma 2.1 exhaustively: all integer X in [1, 3000], c in grid, L = floor(c sqrt X) (exact isqrt):
   max_{d != 0} r_I(d) <= 1 + floor(c^2).
R3 Prop 2.3(2): all n <= 2*10^5: #{(x,y) >= 0 : x^2+y^2 = n, |x-y| <= c n^(1/4)} <= 2(c^2/sqrt2 + 1).
R4 Cor 2.2: all q <= 30, K <= 40, a in [ceil((qK/2)^(4/3)), that + 3000]: #squares <= sqrt(6K).
R5 Lemma 4.2(b) 'Consequently' clause WITHOUT the hidden hypothesis psi <= pi/4 (small R):
   all N <= 3000, w in a list, kappa, C in grids: count of points within (kappa + C/2)/sqrt R of arg w
   versus 2|w|(kappa+C)^2 + 2.  (Floating point only in the angle test; bound is far from tight.)
R6 at alpha = 1/4 exactly two cluster points share a residue class: z = u^2, w = i u ubar,
   u = A + (A-1) i: z = w = 0 mod u, Norm(u) = Norm(N)^(1/4) with N = Norm(u)^2, |z-w| = sqrt2 N^(1/4).
R7 the family u = A + i (claimed: attains the bound up to sqrt2) needs A even.
"""
import math, itertools
from fractions import Fraction

# R2
worst = 0
for X in range(1, 3001):
    for c2num in [1, 2, 4, 9, 16]:          # c^2 in {1/4, 1/2, 1, 9/4, 4}
        c2 = Fraction(c2num, 4)
        L = math.isqrt(int(c2 * X)) if c2 * X == int(c2 * X) else math.isqrt(math.floor(c2 * X))
        # largest L with L^2 <= c^2 X
        while (L+1)**2 <= c2 * X: L += 1
        while L*L > c2 * X: L -= 1
        xs = range(X, X + L + 1)
        cnt = {}
        for x in xs:
            for xp in xs:
                if x > xp:
                    d = x*x - xp*xp
                    cnt[d] = cnt.get(d, 0) + 1
        mx = max(cnt.values()) if cnt else 0
        bound = 1 + math.floor(c2)
        assert mx <= bound, (X, c2, mx)
        worst = max(worst, mx - bound)
print(f"R2 Lemma 2.1: all X <= 3000, c^2 in {{1/4,1/2,1,9/4,4}}: no violation (max excess {worst})")

# R3
viol = 0; tight = 0
for n in range(1, 200001):
    pts = []
    r = math.isqrt(n)
    for x in range(0, r + 1):
        y2 = n - x*x; y = math.isqrt(y2)
        if y*y == y2: pts.append((x, y))
    for c in [0.5, 1.0, 2.0]:
        k = sum(1 for (x, y) in pts if (x - y)**4 <= c**4 * n)
        b = 2*(c*c/math.sqrt(2) + 1)
        if k > b: viol += 1
        if k > 2*(c*c/math.sqrt(2)): tight += 1
print(f"R3 Prop 2.3(2) (with the sharper constant 2(c^2/sqrt2+1)): n <= 2e5, c in {{.5,1,2}}: violations {viol}")

# R4
worst = 0.0
for q in range(1, 31):
    for K in range(1, 41):
        a0 = math.ceil((Fraction(q*K, 2))**Fraction(4, 3)) if False else None
        # exact ceiling of (qK/2)^(4/3): smallest a with a^3 >= (qK/2)^4
        t = Fraction(q*K, 2)**4
        a0 = max(1, int(round(float(t) ** (1/3))) - 2)
        while Fraction(a0)**3 < t: a0 += 1
        for a in range(a0, a0 + 3000):
            lo = math.isqrt(a - 1) + 1 if a > 0 else 0
            hi = math.isqrt(a + q*(K-1))
            c = sum(1 for y in range(lo, hi + 1) if (y*y - a) % q == 0)
            assert c*c <= 6*K, (a, q, K, c)
            worst = max(worst, c / math.sqrt(6*K))
print(f"R4 Cor 2.2: q<=30, K<=40, 3000 values of a from the threshold: no violation; worst ratio {worst:.3f}")

# R5
ws = [(1, 0), (1, 1), (2, 1), (1, 2), (3, 1)]
viol = 0; big_psi = 0; tested = 0
for N in range(1, 3001):
    r = math.isqrt(N); pts = []
    for x in range(-r, r + 1):
        y2 = N - x*x; y = math.isqrt(y2)
        if y*y == y2:
            pts.append((x, y))
            if y: pts.append((x, -y))
    if not pts: continue
    R = math.sqrt(N)
    for w in ws:
        aw = math.atan2(w[1], w[0])
        for kappa in [0.0, 0.5, 1, 2, 4, 8]:
            for C in [0.1, 0.5, 1, 2, 4]:
                psi = (kappa + C/2)/math.sqrt(R)
                near = sum(1 for z in pts
                           if abs(((math.atan2(z[1], z[0]) - aw + math.pi) % (2*math.pi)) - math.pi) <= psi + 1e-12)
                b = 2*math.sqrt(gn := w[0]**2 + w[1]**2)*(kappa + C)**2 + 2
                tested += 1
                if psi > math.pi/4: big_psi += 1
                if near > b: viol += 1
print(f"R5 Lemma 4.2(b) consequence: {tested} cases ({big_psi} with psi > pi/4, outside the proof): violations {viol}")

# R6
def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
for A in [2, 3, 10, 101, 1000]:
    u = (A, A - 1); n = u[0]**2 + u[1]**2
    z = gmul(u, u); w = gmul((0, 1), gmul(u, (u[0], -u[1])))
    N = n*n
    assert z[0]**2 + z[1]**2 == N and w[0]**2 + w[1]**2 == N
    d = (z[0] - w[0], z[1] - w[1])
    num = gmul(d, (u[0], -u[1]))
    assert num[0] % n == 0 and num[1] % n == 0           # u | z - w (indeed u | z and u | w)
    assert d[0]**2 + d[1]**2 == 2*n                      # |z - w|^2 = 2|u|^2 = 2 N^(1/2)
print("R6 alpha = 1/4: z=u^2, w=i*u*ubar share the class 0 mod u, Norm(u)=Norm(N)^(1/4), |z-w|=sqrt2 N^(1/4): verified A=2,3,10,101,1000")

# R7
bad = [A for A in range(1, 50) if (A*A + 1) % 2 == 0]
print(f"R7 u=A+i: Norm(u)=A^2+1 is even (so (1+i) | gcd(u,ubar)) for odd A, e.g. {bad[:5]}; the family needs A even")
