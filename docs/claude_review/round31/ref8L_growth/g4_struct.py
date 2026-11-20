#!/usr/bin/env python3
"""Structured prime sets for F_delta(P):
 (a) P = union of prime factors of t^2+1, t in {t0, t0+2, ..., t0+2(L-1)} (t even): the b_t = t+i
     have |eps b - conj b|^2 = 4 (y = 1), so each b_t is a candidate solution;
 (b) P = union of prime factors of the norms of several unrelated two-mark critical clusters.
Exact counts (meet in the middle + exact verification) up to h <= B with log B = 1.1 log prod P + 10.
"""
import math, sys
sys.path.insert(0, '.')
from fractions import Fraction
from sympy import factorint
from g2_F import count
from g3_heur import expected
def run(ps, label, deltas=('49/100', '3/10'), gammas=(0.25, 0.4, 0.49)):
    ps = sorted(set(ps))
    lp = sum(math.log(p) for p in ps)
    logB = lp + 12
    for d in deltas:
        res, sols, nex, n1, n2 = count(ps, logB, Fraction(d), list(gammas), verify_all=True)
        # how many of the solutions are single b_t (y=+-1 in some sector) / maximal height
        mx = max((lh for e, lh in sols), default=0)
        mod = {g: round(expected(ps, float(Fraction(d)), g, logB, dt=0.002), 1) for g in gammas}
        tail = round(expected(ps, float(Fraction(d)), gammas[0], 3000.0, dt=0.01, logY0=logB), 2)
        print('%s s=%d log prodP=%.1f delta=%s F(by gamma)=%s model=%s model_tail_beyond_B=%s max log h/log prodP=%.3f' % (label, len(ps), lp, d, res, mod, tail, mx/lp))
mode = sys.argv[1]
if mode == 'a':
    for t0, Ls in ((1000, (1, 2, 3, 4, 6)), (30000, (1, 2, 3, 4)), (10**6, (1, 2, 3, 4, 5))):
        for L in Ls:
            ps = []
            for j in range(L):
                t = t0 + 2*j
                ps += list(factorint(t*t+1))
            run(ps, '(a) t0=%d L=%d' % (t0, L))

if mode == 'b':
    # unrelated two-mark clusters: N_j = n_j (a_j^2+1), marks g(a+i), g(a-i) with |g|^2 = n_j ~ a^2/4 (C=1);
    # P = union of primes of all N_j.  Use a_j even with a_j^2+1 having only primes = 1 mod 4 (automatic),
    # n_j a split prime near a_j^2/4.
    from sympy import nextprime
    import random
    random.seed(3)
    for r in (1, 2, 3, 4, 5):
        ps = []
        for j in range(r):
            a = 2*random.randint(500, 3000)
            q = nextprime(a*a//4)
            while q % 4 != 1: q = nextprime(q)
            ps += list(factorint(a*a+1)) + [q]
        run(ps, '(b) r=%d clusters' % r)
if mode == 'c':
    # P = {t^2+1 prime}, t even near 10^6: no small primes, so the random-model count is ~0 and
    # every solution is "structural" (products of b_t^{+-1}).
    from sympy import isprime
    ts = []
    t = 10**6
    while len(ts) < 14:
        if isprime(t*t+1): ts.append(t)
        t += 2
    for L in (2, 4, 6, 8, 10, 12, 14):
        run([t*t+1 for t in ts[:L]], '(c) L=%d t in [%d,%d]' % (L, ts[0], ts[L-1]), gammas=(0.1, 0.25, 0.4))
