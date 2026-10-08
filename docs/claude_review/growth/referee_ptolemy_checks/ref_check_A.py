"""Referee check of Theorem A (general Gaussian Ptolemy with matching coefficients), on ACTUAL
circles x^2+y^2=N with N a product of split prime powers, all four units allowed (mixed units).
Quadruples are taken in true angular (cyclic) order, including consecutive windows (short arcs).
Checks (exact):
 (1) z_i - z_j = g_ij v_ij,  v_ij = eps_i h_ij - eps_j conj(h_ij)
 (2) g12 g34 = Gamma X1, g13 g24 = Gamma X2, g14 g23 = Gamma X3   in Z[i]
 (3) X2 v13 v24 = X1 v12 v34 + X3 v14 v23                        in Z[i]
 (4) the three v-products are positive integer multiples of one primitive omega (cyclic order)
 (5) |v_ij|^2 formulas by unit ratio
"""
import sys, math, itertools, random
from ref_common import *

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 11)

def build_points(primes, exps):
    pis = {p: two_squares(p) for p in primes}
    pts = []
    for alloc in itertools.product(*[range(e+1) for e in exps]):
        z = (1, 0)
        for p, e, a in zip(primes, exps, alloc):
            z = gmul(z, gmul(gpow(pis[p], a), gpow(gconj(pis[p]), e - a)))
        for u, eps in enumerate(UNITS):
            pts.append((gmul(eps, z), u, alloc))
    return pis, pts

def gcd_g(a):
    # content of Gaussian integer as rational gcd of coords (for primitivity of omega)
    return math.gcd(a[0], a[1])

stats = dict(circles=0, quads=0, mixed=0, windows=0, nonunitX=0)
for trial in range(60):
    r = random.randint(1, 4)
    primes = random.sample(SPLIT[:12], r)
    exps = [random.randint(1, 3) for _ in primes]
    N = math.prod(p**e for p, e in zip(primes, exps))
    pis, pts = build_points(primes, exps)
    pts.sort(key=lambda t: math.atan2(t[0][1], t[0][0]))
    n = len(pts)
    if n < 4:
        continue
    stats['circles'] += 1
    cand = []
    for s in range(n):  # consecutive windows (genuinely short arcs)
        cand.append(tuple((s + j) % n for j in range(4)))
    for _ in range(150):
        cand.append(tuple(sorted(random.sample(range(n), 4))))
    for idx in cand:
        Q = [pts[i] for i in idx]
        if len(set(q[0] for q in Q)) < 4:
            continue
        zs = [q[0] for q in Q]; eps = [UNITS[q[1]] for q in Q]; al = [q[2] for q in Q]
        def gh(i, j):
            g = (1, 0); h = (1, 0)
            for k, p in enumerate(primes):
                pi = pis[p]; e = exps[k]; ai, aj = al[i][k], al[j][k]
                g = gmul(g, gmul(gpow(pi, min(ai, aj)), gpow(gconj(pi), e - max(ai, aj))))
                if ai > aj:
                    h = gmul(h, gpow(pi, ai - aj))
                elif aj > ai:
                    h = gmul(h, gpow(gconj(pi), aj - ai))
            return g, h
        G = {}; V = {}
        for i, j in itertools.combinations(range(4), 2):
            g, h = gh(i, j)
            v = gsub(gmul(eps[i], h), gmul(eps[j], gconj(h)))
            assert gmul(eps[i], gmul(g, h)) == zs[i] and gmul(eps[j], gmul(g, gconj(h))) == zs[j]
            assert gsub(zs[i], zs[j]) == gmul(g, v)                      # (1)
            G[i, j] = g; V[i, j] = v
            # (5)
            ratio = gdivexact(eps[j], eps[i])
            x, y = h
            if ratio == (1, 0):
                assert gnorm(v) == 4*y*y
            elif ratio == (-1, 0):
                assert gnorm(v) == 4*x*x
            elif ratio == (0, 1):
                assert gnorm(v) == 2*(x - y)**2
            else:
                assert gnorm(v) == 2*(x + y)**2
        # Gamma and X_k
        Gam = (1, 0); X = [1, 1, 1]
        matchings = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
        for k, p in enumerate(primes):
            pi = pis[p]; e = exps[k]
            vals = [al[i][k] for i in range(4)]
            b = sorted(vals)
            Gam = gmul(Gam, gmul(gpow(pi, b[0] + b[1]), gpow(gconj(pi), 2*e - b[2] - b[3])))
            mp = b[2] - b[1]
            if mp > 0:
                low = sorted(range(4), key=lambda i: vals[i])[:2]
                low = tuple(sorted(low))
                for t, (A, B) in enumerate(matchings):
                    if A == low or B == low:
                        X[t] *= p**mp
        for t, (A, B) in enumerate(matchings):
            assert gmul(G[A], G[B]) == gscal(X[t], Gam), "common factor (2)"      # (2)
        lhs = gscal(X[1], gmul(V[0, 2], V[1, 3]))
        rhs = gadd(gscal(X[0], gmul(V[0, 1], V[2, 3])), gscal(X[2], gmul(V[0, 3], V[1, 2])))
        assert lhs == rhs, "Gaussian Ptolemy (3)"
        # (4): all three products on one ray
        prods = [gmul(V[0, 1], V[2, 3]), gmul(V[0, 2], V[1, 3]), gmul(V[0, 3], V[1, 2])]
        c = gcd_g(prods[1]); om = (prods[1][0]//c, prods[1][1]//c)
        for pr in prods:
            cc = gcd_g(pr)
            assert (pr[0]//cc, pr[1]//cc) == om, "same primitive ray (4)"
        assert math.gcd(X[0], X[1]) == math.gcd(X[0], X[2]) == math.gcd(X[1], X[2]) == 1
        stats['quads'] += 1
        stats['mixed'] += len(set(q[1] for q in Q)) > 1
        stats['nonunitX'] += max(X) > 1
    stats['windows'] += 1
print("Theorem A independent check passed:", stats)
