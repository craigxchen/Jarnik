"""Lemma 3.1 (cofactor) exact check and chain deficits (Prop 3.5) on genuine tuples and on the split Pell points.
(1) random genuine 5-tuples (with multi-layer primes): for random supports S in M and random f supported on S,
    D_S(zeta) | F(zeta) in Z[i] and N(F)/N(D) = N(f(y)) * H_S(y)^2 with log H_S = 1/2 sum_p width_S(a(p)) log p,
    checked as the exact integer identity N(F) * prod_p p^{width_S(a(p))} = N(c)... i.e. N(c) = N(F)/N(D).
(2) chain deficit Delta_S on the same tuples; it vanishes on simple profiles.
(3) split Pell points: chain deficits for S = {0, m} with m = e_a - e_b + e_c - e_d (all 4-tuples of indices) are O(1).
"""
import random
from math import log
from itertools import combinations, permutations
from gauss_lib import *
from cox_form import build_tuple, random_split_primes, CUTS

def width(S, v):
    vals = [sum(m[j]*v[j] for j in range(5)) for m in S]
    return max(vals) - min(vals)

def check_lemma(rng, trials=200):
    ok = 0; tot = 0; maxdef = 0
    for _ in range(trials):
        primes = random_split_primes(rng.randint(2, 5), 5, 120, rng)
        eps, prof, zeta = build_tuple(primes, rng, emax=4)
        if len(set(zeta)) < 5:
            continue
        # random support: 2..4 points in M with small entries
        S = set()
        while len(S) < rng.randint(2, 4):
            v = [rng.randint(-1, 1) for _ in range(4)]
            S.add(tuple(v + [-sum(v)]))
        S = sorted(S)
        coef = {m: complex(rng.randint(-3, 3), rng.randint(-3, 3)) for m in S}
        coef = {m: (int(c.real), int(c.imag)) for m, c in coef.items()}
        if all(c == (0, 0) for c in coef.values()):
            continue
        m0 = [min(m[j] for m in S) for j in range(5)]
        F = (0, 0)
        for m in S:
            term = coef[m]
            for j in range(5):
                term = gmul(term, gpow(zeta[j], m[j] - m0[j]))
            F = (F[0] + term[0], F[1] + term[1])
        D = (1, 0)
        logH2 = 0.0; Nwidth = 1; simplecut = 0.0
        for p in primes:
            pi = gauss_prime(p); e, a = prof[p]
            mn = min(sum(m[j]*a[j] for j in range(5)) for m in S)
            mx = max(sum(m[j]*a[j] for j in range(5)) for m in S)
            ea = [e - x for x in a]
            D = gmul(D, gpow(pi, mn - sum(m0[j]*a[j] for j in range(5))))
            D = gmul(D, gpow(gconj(pi), -mx - sum(m0[j]*ea[j] for j in range(5))))
            Nwidth *= p**(mx - mn)
            # simple (cut-sum) value
            for t in range(1, e+1):
                v = [1 if a[j] >= t else 0 for j in range(5)]
                simplecut += width(S, v) * log(p)
            logH2 += (mx - mn) * log(p)
        tot += 1
        if F == (0, 0):
            ok += 1; continue
        if not gdivides(D, F):
            print('DIVISIBILITY FAIL', primes); continue
        c = gexactdiv(F, D)
        # |c|^2 = |f(y)|^2 H_S^2, |F| = |f(y)| R''^{-sum m0}:  N(c) = N(F)/N(D), N(D) = N'^{-sum m0} / prod p^{width}
        Nprime = gnorm(zeta[0])
        lhs = gnorm(c) * (Nprime ** (-sum(m0)))
        rhs = gnorm(F) * Nwidth
        ok += (lhs == rhs)
        maxdef = max(maxdef, simplecut - logH2)
    print('(1) Lemma 3.1: %d/%d exact identities hold; max chain deficit seen = %.2f' % (ok, tot, maxdef))

def pell_deficits():
    from split_pell import F as PF, split_cox
    from sympy import factorint
    from math import gcd
    x, z = 1, 11
    pts = []
    for _ in range(7):
        x, z = 15*x + 2*z, 15*z + 112*x
        pts.append((x, z))
    for (x, z) in pts:
        u = PF(x, z)
        g = 0
        for y in u: g = gcd(g, y)
        u = [y//g for y in u]
        fac = {}
        for y in u:
            fac.update(factorint(y))
        worst = 0.0
        for idx in permutations(range(5), 4):
            m = [0]*5; m[idx[0]] += 1; m[idx[1]] -= 1; m[idx[2]] += 1; m[idx[3]] -= 1
            S = [tuple([0]*5), tuple(m)]
            dsum = 0.0
            for p in fac:
                v = [0]*5
                for j in range(5):
                    y = u[j]
                    while y % p == 0:
                        y //= p; v[j] += 1
                full = width(S, v)
                cut = sum(width(S, [1 if v[j] >= t else 0 for j in range(5)]) for t in range(1, max(v)+1))
                dsum += (cut - full) * log(p)
            worst = max(worst, dsum)
        print('(3) Pell x=%-12d max chain deficit over S={0,e_a-e_b+e_c-e_d}: %.2f   (w ~ log x/2 = %.1f)' % (x, worst, log(x)/2))

if __name__ == '__main__':
    check_lemma(random.Random(31), 400)
    pell_deficits()
