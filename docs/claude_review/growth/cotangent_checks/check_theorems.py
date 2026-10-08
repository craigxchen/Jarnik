"""Exact checks of the master identity (S2), Theorem A's ingredients and
Theorem B (five-point rigidity) on actual circles.

(S2)  N^Q * prod_{i<j}(X_i-X_j)^2 * prod_e g_e^2 eps_e = Delta * prod_i (X_i^2+L^2)^k
      Q = floor((k+1)^2/4), Delta = prod_p p^{sum_tau (S_tau-M/2)^2 * 4}/... (see text)
      -- we check the equivalent integer form with 4*Delta exponent handled exactly.
(A)   level capacity: at split p, #points per allocation level <= (p-1) p^{v_p(L0)}.
(B)   for M=5: sigma * prod s_ij <= K^10/128 with K = maxchord / N^(1/5);
      6 | L0 for every 5-point cluster.
"""
import math, random
from fractions import Fraction
from itertools import combinations
from cot_lib import (gmul, gconj, gnorm, gsub, gexact, ggcd, ggcd_all, gdivides,
                     lcm_all, factor, vp, split_prime_pi, edge_norm, pair_cot,
                     all_edge_lcm, circle_points, angle, half_cot)
from check_dictionary import chart, SPLIT


def alloc(z, pi):
    t = 0
    while gdivides(pi, z):
        z = gexact(z, pi); t += 1
    return t


def check(cl, N, st):
    g = ggcd_all(cl)
    cl = [gexact(z, g) for z in cl]
    N = N // gnorm(g)
    L, X, cots = chart(cl, N)
    M = len(cl); k = M - 1
    Q = (M * M) // 4
    # cut deficit: 4*D_p = sum_tau (2 S_tau - M)^2 (integer)
    fourD = {}
    levels = {}
    for p, e in factor(N).items():
        pi = split_prime_pi(p)
        al = [alloc(z, pi) for z in cl]
        assert min(al) == 0 and max(al) == e
        fourD[p] = sum((2 * sum(1 for a in al if a < tau) - M) ** 2 for tau in range(1, e + 1))
        levels[p] = al
    # (S2) exact identity: prod_e n_e = N^Q / Delta, Delta = prod p^{D_p}; use 4D to stay integral
    prod_n = 1
    for i, j in combinations(range(M), 2):
        c = cots[(i, j)]
        a, s = c.numerator, c.denominator
        eps = 2 if (a % 2 and s % 2) else 1
        prod_n *= (a * a + s * s) // eps
    # exact cut identity in integral form: (prod_e n_e)^4 * prod_p p^{4 D_p} = N^{M^2}
    den4 = 1
    for p in fourD:
        den4 *= p ** fourD[p]
    assert prod_n ** 4 * den4 == N ** (M * M)
    # cotangent form: prod_e n_e = prod_i (X_i^2+L^2)^k / (prod (X_i-X_j)^2 prod g_e^2 eps_e)
    num = 1
    for x in X:
        num *= (x * x + L * L) ** k
    den = 1
    for i, j in combinations(range(k), 2):
        den *= (X[i] - X[j]) ** 2
    edges = list(X) + [pair_cot(X[i], X[j], L) for i, j in combinations(range(k), 2)]
    for q in edges:
        gg = math.gcd(q, L)
        a, b = q // gg, L // gg
        eps = 2 if (a % 2 and b % 2) else 1
        den *= gg * gg * eps
    assert num % den == 0 and num // den == prod_n
    # (A) level capacity
    for p, al in levels.items():
        cap = (p - 1) * p ** vp(L, p)
        for lev in set(al):
            assert al.count(lev) <= cap
    st['n'] += 1
    # (B) five points
    if M == 5:
        assert L % 6 == 0, L
        sigma = 1
        for p, e in factor(N).items():
            al = levels[p]
            for tau in range(1, e + 1):
                lo = sum(1 for a in al if a < tau)
                if lo in (1, 4):
                    sigma *= p
        prod_s = 1
        for c in cots.values():
            prod_s *= c.denominator
        maxch = max(math.sqrt(gnorm(gsub(cl[i], cl[j]))) for i, j in combinations(range(5), 2))
        K = maxch / N ** 0.2
        assert sigma * prod_s <= K ** 10 / 128 * (1 + 1e-9), (sigma, prod_s, K)
        st['five'] += 1
        st['minK'] = min(st['minK'], K)


def main(seed=3, trials=1500):
    random.seed(seed)
    st = {'n': 0, 'five': 0, 'minK': 1e9}
    for t in range(trials):
        r = random.randint(1, 5)
        ps = random.sample(SPLIT, r)
        pp = [(p, random.randint(1, 3)) for p in ps]
        N = 1
        for p, e in pp:
            N *= p ** e
        if N > 10 ** 14:
            continue
        pts = sorted(circle_points(pp), key=angle)
        n = len(pts)
        for _ in range(4):
            k = random.choice([2, 3, 4, 4, 4, 5, 6])
            if k + 1 > n // 4:
                continue
            s = random.randrange(n)
            cl = [pts[(s + u) % n] for u in range(k + 1)]
            span = (angle(cl[-1]) - angle(cl[0])) % (2 * math.pi)
            if span >= math.pi / 2 or span == 0:
                continue
            check(cl, N, st)
    print('clusters: %d, five-point clusters: %d, min K (maxchord/N^(1/5)) seen: %.4f'
          % (st['n'], st['five'], st['minK']))
    print('ALL THEOREM CHECKS PASSED')


if __name__ == '__main__':
    main()
