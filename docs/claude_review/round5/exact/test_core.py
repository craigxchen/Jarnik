"""Unit tests of core.py against brute force (exact).

(1) |G_T| from the diagonalisation equals the brute-force count of mu in (1/D)Z^{M-1} cap [0,1)^{M-1}
    with Q mu integral (D = lcm of the d_i, checked also with D' = 2D to make sure nothing is missed).
(2) every element produced by Group.element lies in L_T, distinct n give distinct classes.
(3) kappa_of(eta) realises the homomorphism: mu.kappa = h_eta(mu) mod 1 for all elements.
(4) the exact left inverse satisfies A Q = I and the coercivity bounds hold on random lambda.
(5) verify() agrees with a naive enumeration on small cases.
"""
import random, itertools, math
from fractions import Fraction
import numpy as np
from core import (Group, Qmat, all_patterns, frac_part, exact_left_inverse, left_inverse_lp, coercivity,
                  lcm_list)
from verify import verify

random.seed(7)
ok = True
ntests = 0
for trial in range(400):
    M = random.choice([3, 4, 5, 6])
    pats = all_patterns(M)
    t = random.randint(M - 1, min(len(pats), M + 3))
    T = random.sample(pats, t)
    G = Group(M, T)
    Q = Qmat(M, T)
    n = M - 1
    if not G.full_rank:
        assert np.linalg.matrix_rank(Q.astype(float)) < n
        continue
    ntests += 1
    D = lcm_list(G.invs) if G.invs else 1
    for DD in (D, 2 * D):
        if DD ** n > 300000:
            continue
        cnt = 0
        for v in itertools.product(range(DD), repeat=n):
            if all(sum(Q[i, j] * v[j] for j in range(n)) % DD == 0 for i in range(len(T))):
                cnt += 1
        if cnt != G.order:
            print('ORDER MISMATCH', M, T, G.invs, cnt, DD); ok = False
    seen = set()
    for nv in G.elements():
        mu = G.element(nv)
        for i in range(len(T)):
            s = sum(mu[j] for j in range(n) if Q[i, j])
            if s.denominator != 1:
                print('NOT IN L_T'); ok = False
        key = tuple(frac_part(mu))
        if key in seen:
            print('DUPLICATE CLASS'); ok = False
        seen.add(key)
    for _ in range(3):
        eta = [random.randrange(d) for d in G.invs]
        kap = G.kappa_of(eta)
        assert G.eta_of_kappa(kap) == [e % d for e, d in zip(eta, G.invs)]
        for nv in G.elements():
            mu = G.element(nv)
            if sum(m * k for m, k in zip(mu, kap)) % 1 != G.hom_value(eta, nv):
                print('HOM MISMATCH'); ok = False
    om = [Fraction(random.randint(4, 12), 4) for _ in T]
    A = exact_left_inverse(Q, left_inverse_lp(Q, [float(w) for w in om]))
    k_inf, k_1 = coercivity(Q, om, A)
    for _ in range(50):
        mu = [Fraction(random.randint(-20, 20), random.randint(1, 6)) for _ in range(n)]
        lam = [-sum(mu)] + mu
        N = sum(w * abs(sum(mu[j] for j in range(n) if Q[i, j])) for i, w in enumerate(om))
        if max(abs(v) for v in lam) > k_inf * N or sum(abs(v) for v in lam) > k_1 * N:
            print('COERCIVITY FAIL'); ok = False
print('group/hom/coercivity tests on', ntests, 'full-rank supports:', 'PASS' if ok else 'FAIL')

# (5) verify vs naive enumeration of lambda with small entries
ok5 = True
nt5 = 0
for trial in range(60):
    M = random.choice([4, 5])
    pats = all_patterns(M)
    T = random.sample(pats, random.randint(M - 1, M + 2))
    G = Group(M, T)
    if not G.full_rank:
        continue
    om = [Fraction(random.randint(4, 8), 4) for _ in T]
    eta = [random.randrange(d) for d in G.invs]
    kap = G.kappa_of(eta)
    rho0 = Fraction(random.randint(0, 2), 2)
    rho1 = Fraction(random.randint(0, 1), 4)
    r = verify(M, T, om, kap, rho0, rho1)
    if r['ok'] is None:
        continue
    nt5 += 1
    # naive: lambda with denominators dividing lcm(invs), entries in [-3,3]
    D = lcm_list(G.invs) if G.invs else 1
    Q = Qmat(M, T)
    n = M - 1
    best = None
    rng = range(-3 * D, 3 * D + 1)
    for v in itertools.product(rng, repeat=n):
        if not any(v):
            continue
        mu = [Fraction(a, D) for a in v]
        lam_vals = [sum(mu[j] for j in range(n) if Q[i, j]) for i in range(len(T))]
        if any(x.denominator != 1 for x in lam_vals):
            continue
        if sum(m * k for m, k in zip(mu, kap)) % 1 != 0:
            continue
        lam = [-sum(mu)] + mu
        if max(abs(x) for x in lam) > 3:
            continue
        s = sum(w * abs(x) for w, x in zip(om, lam_vals)) - sum(om) / 2 - rho0 - rho1 * sum(abs(x) for x in lam) / 2
        if best is None or s < best:
            best = s
    if best is not None and best < r['min_slack']:
        print('VERIFY MISSED', M, T, om, kap, rho0, rho1, best, r['min_slack']); ok5 = False
    if best is not None and r['min_slack'] < best and r['Binf'] <= 3:
        print('VERIFY FOUND SMALLER INSIDE NAIVE BOX?', best, r['min_slack']); ok5 = False
print('verify vs naive on', nt5, 'cases:', 'PASS' if ok5 else 'FAIL')
