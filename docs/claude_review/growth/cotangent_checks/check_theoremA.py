"""Sanity check of Theorem A's inequality on actual clusters:
 sum_{split p} log p * (M-4c_p)_+^3/(12 c_p) <= W M/4 + binom(M,2) log(C^2/2),
 c_p=(p-1) p^{v_p(L0)}, C = max chord / N^(1/4), plus the intermediate
 inequality sum_p D_p log p <= W M/4 + binom(M,2) log(C^2/2)."""
import math, random
from itertools import combinations
from cot_lib import (gnorm, gsub, gexact, ggcd_all, gdivides, factor, vp, split_prime_pi,
                     circle_points, angle)
from check_dictionary import chart, SPLIT

random.seed(7)
n_ok = 0
worst = -1e9
for _ in range(3000):
    r = random.randint(1, 5)
    pp = [(p, random.randint(1, 4)) for p in random.sample([q for q in SPLIT if q < 60], r)]
    N = 1
    for p, e in pp:
        N *= p ** e
    if N > 10 ** 13:
        continue
    pts = sorted(circle_points(pp), key=angle)
    n = len(pts)
    if n // 4 < 3:
        continue
    M = random.randint(3, min(12, n // 4))
    s = random.randrange(n)
    cl = [pts[(s + t) % n] for t in range(M)]
    span = (angle(cl[-1]) - angle(cl[0])) % (2 * math.pi)
    if span >= math.pi / 2:
        continue
    g = ggcd_all(cl)
    cl = [gexact(z, g) for z in cl]
    Np = gnorm(cl[0])
    if Np == 1:
        continue
    L, X, _ = chart(cl, Np)
    W = math.log(Np)
    C = max(math.sqrt(gnorm(gsub(a, b))) for a, b in combinations(cl, 2)) / Np ** 0.25
    rhs = W * M / 4 + M * (M - 1) / 2 * math.log(C * C / 2)
    lhsD = 0.0
    lhsA = 0.0
    for p, e in factor(Np).items():
        pi = split_prime_pi(p)
        al = []
        for z in cl:
            t = 0
            while gdivides(pi, z):
                z = gexact(z, pi); t += 1
            al.append(t)
        D = sum((sum(1 for a in al if a < tau) - M / 2) ** 2 for tau in range(1, e + 1))
        lhsD += D * math.log(p)
        c = (p - 1) * p ** vp(L, p)
        lhsA += math.log(p) * max(M - 4 * c, 0) ** 3 / (12 * c)
    assert lhsA <= lhsD + 1e-9
    assert lhsD <= rhs + 1e-6, (lhsD, rhs)
    worst = max(worst, lhsD - rhs)
    n_ok += 1
print('Theorem A chain verified on %d clusters (max of LHS_D - RHS = %.3f)' % (n_ok, worst))
