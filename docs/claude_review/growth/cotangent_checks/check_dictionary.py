"""Exact check of the integer-cotangent dictionary on actual circles.

For random primitive Gaussian circles x^2+y^2=N (N a product of split prime
powers) we take k+1 angularly consecutive lattice points, divide by their
Gaussian gcd, and verify every identity of Section 1 of cotangent.md:

  (D1) chart: X_j = L_0 cot(theta_j/2) integral, Q_ij = L_0 cot(theta_ij/2)
       equal to (X_i X_j + L^2)/(X_i - X_j);
  (D2) N_prim = all-edge lcm of reduced edge norms;
  (D3) the primitive realisation reproduces the relative phases;
  (D4) angular span = 2 arctan(L/A);
  (D5) chord formula |z_i-z_j|^2 eps n = 4 N s^2, n_ij = N/Norm(gcd);
  (D6) good-prime valuation dictionary (p not dividing 2L);
  (D7) pair gcd formula Norm(gcd(H_i,H_j)) = |X_i-X_j| d;
  (D8) endpoint swap preserves N, L and A;
  (D9) cut identity: prod_e n_e = prod_p p^{sum_tau S_tau (M-S_tau)}.
"""
import math
import random
import sys
from fractions import Fraction
from itertools import combinations

from cot_lib import (gmul, gconj, gnorm, gsub, gexact, ggcd, ggcd_all, gdivides,
                     lcm_all, factor, vp, split_prime_pi, edge_norm, pair_cot,
                     is_clique, all_edge_lcm, primitive_tuple, circle_points,
                     angle, half_cot)

SPLIT = [p for p in range(5, 120) if p % 4 == 1 and all(p % q for q in range(2, p))]


def chart(cluster, N):
    """cluster: list of Gaussian integers sorted by angle in a short arc.
    Anchor = cluster[0].  Returns L0, X (finite cotangents of cluster[1:])."""
    M = len(cluster)
    cots = {}
    for i, j in combinations(range(M), 2):
        c = half_cot(cluster[i], cluster[j], N)
        cots[(i, j)] = c
    L0 = lcm_all(c.denominator for c in cots.values())
    X = [L0 * cots[(0, j)] for j in range(1, M)]
    assert all(x.denominator == 1 for x in X)
    X = [int(x) for x in X]
    # pair cotangents: angle from i to j, i<j (counterclockwise)
    for i, j in combinations(range(1, M), 2):
        Q = L0 * cots[(i, j)]
        assert Q.denominator == 1
        # cot((t_j - t_i)/2) = (x_i x_j + 1)/(x_i - x_j)
        assert int(Q) == pair_cot(X[i - 1], X[j - 1], L0), (Q, X, L0)
    return L0, X, cots


def check_cluster(cluster, N, stats):
    M = len(cluster)
    g = ggcd_all(cluster)
    cl = [gexact(z, g) for z in cluster]
    Np = N // gnorm(g)
    assert all(gnorm(z) == Np for z in cl)
    L, X, cots = chart(cl, Np)
    k = M - 1
    A = min(X)
    assert X[-1] == A and all(x > 0 for x in X)  # far endpoint has least X
    assert is_clique(X, L)
    # (D2)
    Nl = all_edge_lcm(X, L)
    assert Nl == Np, (Nl, Np)
    # (D3) primitive realisation = cluster up to a common unit
    rows = primitive_tuple(X, L)
    for j in range(M):
        # rows[j]/rows[0] == cl[j]/cl[0]
        assert gmul(rows[j], cl[0]) == gmul(rows[0], cl[j])
    assert gnorm(rows[0]) == Np
    # (D4) span
    span = (angle(cl[-1]) - angle(cl[0])) % (2 * math.pi)
    assert abs(span - 2 * math.atan(L / A)) < 1e-9 * max(1, span)
    Cstar = 2 * Np ** 0.25 * math.atan(L / A)
    stats['Cstar'].append(Cstar)
    # (D5) chords
    for i, j in combinations(range(M), 2):
        c = cots[(i, j)]
        a, s = c.numerator, c.denominator
        eps = 2 if (a % 2 and s % 2) else 1
        n = (a * a + s * s) // eps
        d2 = gnorm(gsub(cl[i], cl[j]))
        assert d2 * eps * n == 4 * Np * s * s
        # eps = 2 exactly when the x-coordinates have different parity
        assert (eps == 2) == ((cl[i][0] - cl[j][0]) % 2 == 1)
        assert n == Np // gnorm(ggcd(cl[i], cl[j]))
        assert Np % n == 0
    # (D6) valuation dictionary at p not dividing 2L
    for p in factor(Np):
        if (2 * L) % p == 0:
            continue
        pi = split_prime_pi(p)
        pib = gconj(pi)
        def vpi(z, w):
            e = 0
            while z != (0, 0) and gdivides(w, z):
                z = gexact(z, w)
                e += 1
            return e
        al = [vpi((x, L), pi) for x in X]
        ga = [vpi((x, L), pib) for x in X]
        assert all(min(a_, g_) == 0 for a_, g_ in zip(al, ga))
        assert vp(Np, p) == max([0] + al) + max([0] + ga)
        # allocation of the primitive points
        alloc = [vpi(z, pi) for z in cl]
        assert alloc[0] == max([0] + ga)
        for j in range(1, M):
            assert alloc[j] == max([0] + ga) + al[j - 1] - ga[j - 1]
            assert vp(X[j - 1] ** 2 + L * L, p) == al[j - 1] + ga[j - 1]
        for i, j in combinations(range(k), 2):
            assert vp(X[i] - X[j], p) == min(al[i], al[j]) + min(ga[i], ga[j])
        stats['goodprimes'] += 1
    # (D7) pair gcd formula
    for i, j in combinations(range(k), 2):
        x, y = X[i], X[j]
        dlt = x - y
        d = math.gcd(math.gcd((x * x + L * L) // dlt, x), math.gcd(L, dlt))
        assert gnorm(ggcd((x, L), (y, L))) == abs(dlt) * d
        assert L % d == 0
    # (D8) endpoint swap
    S = A * A + L * L
    Xs = [A] + [A + S // (x - A) for x in X if x != A]
    assert all(S % (x - A) == 0 for x in X if x != A)
    assert is_clique(Xs, L) and min(Xs) == A
    assert all_edge_lcm(Xs, L) == Np
    # (D9) cut identity
    prod_n = 1
    for i, j in combinations(range(M), 2):
        prod_n *= Np // gnorm(ggcd(cl[i], cl[j]))
    rhs = 1
    for p, e in factor(Np).items():
        pi = split_prime_pi(p)
        def vpi2(z):
            t = 0
            while gdivides(pi, z):
                z = gexact(z, pi)
                t += 1
            return t
        al = [vpi2(z) for z in cl]
        assert min(al) == 0 and max(al) == e
        cut = sum(sum(1 for a_ in al if a_ < tau) * sum(1 for a_ in al if a_ >= tau)
                  for tau in range(1, e + 1))
        rhs *= p ** cut
    assert prod_n == rhs
    stats['clusters'] += 1
    stats['maxk'] = max(stats['maxk'], k)


def main(seed=1, trials=400):
    random.seed(seed)
    stats = {'clusters': 0, 'goodprimes': 0, 'maxk': 0, 'Cstar': []}
    for t in range(trials):
        r = random.randint(1, 4)
        ps = random.sample(SPLIT, r)
        pp = [(p, random.randint(1, 3)) for p in ps]
        N = 1
        for p, e in pp:
            N *= p ** e
        if N > 10 ** 12:
            continue
        pts = circle_points(pp)
        pts.sort(key=angle)
        n = len(pts)
        for _ in range(3):
            k = random.randint(2, min(7, n // 4))
            s = random.randrange(n)
            cl = [pts[(s + t) % n] for t in range(k + 1)]
            # need a short arc (< pi/2) so that the chart is positive
            span = (angle(cl[-1]) - angle(cl[0])) % (2 * math.pi)
            if span >= math.pi / 2 or span == 0:
                continue
            check_cluster(cl, N, stats)
    print('clusters checked:', stats['clusters'], ' good-prime checks:',
          stats['goodprimes'], ' max finite coords:', stats['maxk'])
    print('min C_* seen: %.4f' % min(stats['Cstar']))
    print('ALL DICTIONARY CHECKS PASSED')


if __name__ == '__main__':
    main(*(int(a) for a in sys.argv[1:]))
