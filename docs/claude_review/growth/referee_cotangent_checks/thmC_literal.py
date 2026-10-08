# Referee: Theorem C as literally stated ("prod z_i = eps * N_bal^(M/2) * Phi with Norm(Phi) <= Delta^2",
# N_bal = product of p over balanced layers) versus the corrected statement
# ("prod z_i = eps * N_bal^(M/2) * m * Phi, m = prod_{unbalanced layers} p^min(h,M-h) a positive integer,
#   Norm(Phi) = prod_{layers} p^|M-2h| <= Delta^2").  Exhaustive over even-M cyclic windows, N <= NMAX.
import sys, math
from math import isqrt
sys.path.insert(0, '.')
from rlib import gmul, gnorm, gdivmod_exact, ggcd_list, gauss_prime_above, vpi, angle, factor

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
pts = {}
for x in range(0, isqrt(NMAX) + 1):
    for y in range(0, isqrt(NMAX - x * x) + 1):
        n = x * x + y * y
        if n: pts.setdefault(n, set()).update({(x, y), (-x, y), (x, -y), (-x, -y)})
lit_fail = 0; corr_ok = 0; tested = 0; first = None
for N in sorted(pts):
    L = sorted(pts[N], key=angle); n = len(L)
    for M in (4, 6):
        if n < M: continue
        for st in range(n):
            W = [L[(st + t) % n] for t in range(M)]
            angs = sorted(angle(z) for z in W)
            gaps = [(angs[(t + 1) % M] - angs[t]) % (2 * math.pi) for t in range(M)]
            if max(gaps) <= math.pi + 1e-9: continue
            g = ggcd_list(W); P = [gdivmod_exact(z, g) for z in W]
            Np = gnorm(P[0])
            prod = (1, 0)
            for z in P: prod = gmul(prod, z)
            Nbal = 1; m = 1; normPhi = 1; Delta = 1
            for p, e in factor(Np).items():
                pi = gauss_prime_above(p); lv = [vpi(z, pi) for z in P]
                for tau in range(1, e + 1):
                    h = sum(1 for v in lv if v >= tau)
                    Delta *= p ** ((M // 2 - h) ** 2)
                    normPhi *= p ** abs(M - 2 * h)
                    if 2 * h == M: Nbal *= p
                    else: m *= p ** min(h, M - h)
            tested += 1
            q = gdivmod_exact(prod, (Nbal ** (M // 2), 0))
            if q is None or gnorm(q) > Delta ** 2:
                lit_fail += 1
                span = max(math.dist(a, b) for a in P for b in P) / Np ** 0.25
                if first is None or (M, Np) < (first[0], first[1]):
                    first = (M, Np, P, prod, Nbal, Delta, gnorm(q) if q else None, m, span)
            q2 = gdivmod_exact(prod, (Nbal ** (M // 2) * m, 0))
            if q2 is not None and gnorm(q2) == normPhi and normPhi <= Delta ** 2: corr_ok += 1
            else: print('CORRECTED FAIL', N, W)
print('even windows tested', tested, 'literal statement fails on', lit_fail, 'corrected statement holds on', corr_ok)
M, Np, P, prod, Nbal, Delta, nq, m, span = first
print('smallest literal counterexample: M=%d N=%d points=%s prod=%s N_bal=%d Delta=%d Norm(prod/N_bal^(M/2))=%s Delta^2=%d rational factor m=%d'
      % (M, Np, P, prod, Nbal, Delta, nq, Delta ** 2, m))
