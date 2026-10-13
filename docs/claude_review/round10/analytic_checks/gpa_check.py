# Evidence for Conjecture GPA(h) (analytic.md 5.1).
# Part A (floating point, evidence): random phi; smallest-norm Gaussian prime varpi (Norm != h) with
#   sin(dist(phi - arg varpi, (pi/2)Z)) < (h Norm varpi)^(-1/2).  Associates are covered by the dist mod pi/2.
# Part B (exact integers): lattice directions beta (primitive, h < Norm beta <= 20h): a Gaussian prime on a line
#   Im(alpha conj beta) = m, |m| < sqrt(Norm beta / h)  (m = 0 allowed iff beta is itself prime).
import numpy as np, math, random
from sympy import isprime, primerange
XMAX = 2_000_000
# Gaussian primes: split ones in the open first quadrant (both conjugates), 1+i, inert q (angle 0, norm q^2)
sieve = np.ones(XMAX + 1, dtype=bool); sieve[:2] = False
for i in range(2, int(XMAX**0.5) + 1):
    if sieve[i]: sieve[i*i::i] = False
pts = []
for a in range(1, int(math.isqrt(XMAX)) + 1):
    bmax = int(math.isqrt(XMAX - a*a))
    for b in range(1, bmax + 1):
        n = a*a + b*b
        if sieve[n] and n % 4 == 1: pts.append((n, math.atan2(b, a)))
pts.append((2, math.pi/4))
for q in range(3, int(math.isqrt(XMAX)) + 1, 4):
    if sieve[q]: pts.append((q*q, 0.0))
pts.sort()
norms = np.array([p[0] for p in pts], dtype=float); angs = np.array([p[1] for p in pts])
print("Gaussian primes (classes mod units) up to norm %d: %d" % (XMAX, len(pts)))

def first_witness(phis, h):
    res = np.full(len(phis), np.inf)
    alive = np.arange(len(phis))
    CH = 20000
    for s in range(0, len(pts), CH):
        if len(alive) == 0: break
        nn = norms[s:s+CH]; aa = angs[s:s+CH]
        ok = nn != h
        d = np.abs(((phis[alive][None, :] - aa[:, None]) + np.pi/4) % (np.pi/2) - np.pi/4)
        cond = (np.sin(d) < 1.0/np.sqrt(h*nn)[:, None]) & ok[:, None]
        hit = cond.any(axis=0)
        if hit.any():
            idx = cond.argmax(axis=0)
            res[alive[hit]] = nn[idx[hit]]
            alive = alive[~hit]
    return res

random.seed(1)
for h in (5, 13, 101, 1009):
    phis = np.array([random.uniform(0, np.pi/2) for _ in range(2000)])
    w = first_witness(phis, h)
    fin = w[np.isfinite(w)]
    print("h=%5d random phi (failures are expected near low-height lattice directions): %d/%d have a witness of norm <= %d; median %.0f, max %.0f, max/h %.1f"
          % (h, len(fin), len(w), XMAX, np.median(fin), fin.max(), fin.max()/h))
    # near low-height lattice directions (phi = arg beta + eps, Norm beta <= h): witness norm ~ 1/(h eps^2)?
    # composite primitive beta (no prime parallel to it): 3+4i = (2+i)^2, 4+7i = (2+i)(3+2i)
    for (ba, bb) in ((3, 4), (4, 7)):
        if ba*ba + bb*bb > h: continue
        out = []
        for eps in (1e-2, 3e-3, 1e-3, 3e-4):
            ph = np.array([math.atan2(bb, ba) + eps])
            out.append((eps, first_witness(ph, h)[0]))
        print("   near arg(%d+%di):" % (ba, bb), ", ".join("eps=%g -> %s" % (e, int(x) if np.isfinite(x) else 'none<=XMAX') for e, x in out))

def gp(x, y):
    if x == 0 or y == 0:
        m = abs(x + y); return m % 4 == 3 and isprime(m)
    return isprime(x*x + y*y)

def egcd(u, v):
    if v == 0: return (u, 1, 0)
    g, s, t = egcd(v, u % v); return (g, t, s - (u // v) * t)

print("Part B: lattice directions")
for h in (5, 13, 101):
    worst = (0, None); nb = 0; fails = 0
    H = 20 * h
    for a in range(0, int(math.isqrt(H)) + 1):
        for b in range(1, int(math.isqrt(H)) + 1):
            n = a*a + b*b
            if not (h < n <= H) or math.gcd(a, b) != 1: continue
            nb += 1
            if gp(a, b): continue                      # m = 0 witness (beta itself prime, Norm != h)
            mmax = math.isqrt(n // h) if (math.isqrt(n // h))**2 * h < n else math.isqrt(n // h) - 1
            best = None
            for m in range(-mmax, mmax + 1):
                if m == 0 or m*m*h >= n: continue
                # alpha = x+iy with Im(alpha conj beta) = a*y - b*x = m
                g, s, t = egcd(a, -b)
                y0, x0 = s*m*g, t*m*g
                for k in range(-60, 61):
                    x, y = x0 + k*a, y0 + k*b
                    if (x, y) != (0, 0) and gp(x, y):
                        nrm = x*x + y*y
                        if best is None or nrm < best: best = nrm
            if best is None: fails += 1
            elif best / n > worst[0]: worst = (best / n, (a, b, best))
    print("h=%4d: %d primitive directions with h < Norm <= 20h; failures %d; worst witness norm/Norm beta %.1f at %s"
          % (h, nb, fails, worst[0], worst[1]))
