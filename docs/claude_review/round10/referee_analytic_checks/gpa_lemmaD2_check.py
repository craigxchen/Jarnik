"""Referee re-implementation (independent of analytic_checks/):
(A) GPA(h) at lattice directions: for h in H and every primitive beta = u+iv (u > 0, v >= 0; enough by the
    symmetries i* and conj of the prime set) with h < Norm beta <= 20h, find a Gaussian prime w with Norm w != h
    and h * Im(w conj beta)^2 < Norm beta (exact integers).  Search: lines Im(w conj beta) = m, |m| < sqrt(Nb/h),
    points w = w_m + t beta, |t| <= 2000.
(B) Lemma D.2 hypothesis: for every split p <= PMAX and pi = a+bi, primes w+- with Im(w+- conj pi) = +-1 and
    Re(w+- conj pi) > 0; also checks Norm w+- != p, and (float illustration) that phi just inside the arc
    (arg pi, arg pi + 2 rho+) violates PCG_p through w+.
usage: python3 gpa_lemmaD2_check.py [PMAX]
"""
import math, sys
from sympy import isprime

def gprime(x, y):
    if x == 0 or y == 0:
        m = abs(x + y); return m % 4 == 3 and isprime(m)
    return isprime(x * x + y * y)

def line_point(u, v, m):
    """some w = x+iy with Im(w conj beta) = y u - x v = m (beta = u+iv primitive)."""
    # solve y u - x v = m via extended gcd of (u, -v)
    def eg(a, b):
        if b == 0: return (a, 1, 0)
        g, s, t = eg(b, a % b); return (g, t, s - (a // b) * t)
    g, s, t = eg(u, v)            # u s + v t = g = 1 (up to sign)
    assert abs(g) == 1
    # y u - x v = m  <=  y = s m g, x = -t m g
    y, x = s * m * g, -t * m * g
    assert y * u - x * v == m
    return x, y

# (A)
for h in (29, 53, 197):
    nd = 0; worst = 0; fails = []
    for u in range(1, math.isqrt(20 * h) + 1):
        for v in range(0, math.isqrt(20 * h) + 1):
            Nb = u*u + v*v
            if not (h < Nb <= 20 * h) or math.gcd(u, v) != 1: continue
            nd += 1
            best = None
            mmax = 0
            while h * (mmax + 1) ** 2 < Nb: mmax += 1
            for m in range(-mmax, mmax + 1):
                x0, y0 = line_point(u, v, m)
                for tt in range(0, 2001):
                    for t in ((tt, -tt) if tt else (0,)):
                        x, y = x0 + t * u, y0 + t * v
                        if (x, y) == (0, 0): continue
                        Nw = x*x + y*y
                        if Nw != h and gprime(x, y) and h * m * m < Nb:
                            if best is None or Nw < best: best = Nw
                            break
                    else:
                        continue
                    break
            if best is None: fails.append((u, v))
            else: worst = max(worst, best / Nb)
    print("(A) h=%d: %d primitive directions with h < Norm <= 20h; no witness: %d %s; max witness Norm/Norm(beta) = %.1f"
          % (h, nd, len(fails), fails[:5], worst))

# (B)
PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
cnt = 0; worst_steps = 0; bad = []; arcfail = 0
for p in range(5, PMAX + 1, 4):
    if not isprime(p): continue
    b = 1
    while True:
        a = math.isqrt(p - b*b)
        if a*a + b*b == p: break
        b += 1
    a, b = max(a, b), min(a, b)
    cnt += 1
    for m in (1, -1):
        x0, y0 = line_point(a, b, m)          # Im(w conj pi) = y a - x b = m
        # move along the line to the first point with Re(w conj pi) = x a + y b > 0
        re0 = x0 * a + y0 * b
        t0 = -(re0 // p)                       # re0 + t0 p in (0, p]
        if re0 + t0 * p <= 0: t0 += 1
        found = None
        for s in range(0, 100000):
            x, y = x0 + (t0 + s) * a, y0 + (t0 + s) * b
            if gprime(x, y): found = (x, y, s); break
        if found is None: bad.append((p, m)); continue
        x, y, s = found
        assert y * a - x * b == m and x * a + y * b > 0 and x*x + y*y != p
        worst_steps = max(worst_steps, s)
        if p < 2000 and m == 1:
            # illustration: phi = arg pi + rho (inside the forbidden arc) violates PCG_p via w
            th = math.atan2(b, a); Nw = x*x + y*y
            rho = math.asin((p * Nw) ** -0.5)
            phi = th + rho
            if not abs(math.sin(phi - math.atan2(y, x))) < (p * Nw) ** -0.5: arcfail += 1
print("(B) split p <= %d: %d primes; lines without witness: %d %s; max steps %d; arc illustration failures %d"
      % (PMAX, cnt, len(bad), bad[:5], worst_steps, arcfail))
