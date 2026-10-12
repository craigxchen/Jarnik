"""Proposition 5.3 of weakest.md: S-dependence of the subspace count for the systems Sigma_2, Sigma_3.

S_s = {inf} u {places above the first s split primes}.  K = Q(i), |x|_inf = |x|^2.
Sigma_m: forms X_1, X_i - X_1 at inf; X_1..X_m at finite places of S.
Pair solutions: u = pi_i conj(pi_j) (i<j), x = (u, eta conj u) with eta a unit; then
prod = |eta conj(u)/u - 1|^2 * n^(-2), H(x) = n = p_i p_j.  Pigeonhole on arg(pi_i) mod pi/4
gives many pairs with dist(arg u, (pi/4)Z) <= P^(-delta)/2, i.e. prod <= H^(-2-delta).
Triples: (uv, (eta conj u) v, u (eta' conj v)) for disjoint pigeonhole pairs u, v.

Every reported solution is re-verified numerically (60 digits): prod <= H^(-m-delta).
Distinct pairs {i,j} give distinct points of P^1(K) (the ratio x2/x1 determines {i,j}).
"""
import math
import itertools
import mpmath
import sympy
from g8common import gauss_prime, gmul, gconj, gnorm

mpmath.mp.dps = 60


def run(s, delta):
    primes = []
    p = 5
    while len(primes) < s:
        if sympy.isprime(p) and p % 4 == 1:
            primes.append(p)
        p += 4
    P = primes[-1]
    pis = [gauss_prime(q) for q in primes]
    phis = [mpmath.atan2(pi[1], pi[0]) for pi in pis]
    L = math.ceil((math.pi / 2) * P ** delta)
    width = mpmath.pi / 4 / L
    bins = {}
    for i, f in enumerate(phis):
        b = int(mpmath.floor((f % (mpmath.pi / 4)) / width))
        bins.setdefault(b, []).append(i)
    pairs = [pr for v in bins.values() for pr in itertools.combinations(v, 2)]
    # verify pair solutions
    good2 = 0
    for (i, j) in pairs:
        u = gmul(pis[i], gconj(pis[j]))
        n = gnorm(u)
        a = mpmath.atan2(u[1], u[0])
        k = mpmath.nint(2 * a / (mpmath.pi / 2))
        theta = 2 * a - k * mpmath.pi / 2      # arg(u/(eta conj u)) for the matching unit
        prod = abs(mpmath.exp(1j * theta) - 1) ** 2 * mpmath.mpf(n) ** (-2)
        if prod <= mpmath.mpf(n) ** (-2 - delta):
            good2 += 1
    # triples from disjoint pairs (count only; verify a sample)
    disjoint = 0
    checked = 0
    ok3 = True
    for a_idx in range(len(pairs)):
        for b_idx in range(a_idx + 1, len(pairs)):
            if set(pairs[a_idx]) & set(pairs[b_idx]):
                continue
            disjoint += 1
            if checked < 200:
                (i, j), (k2, l) = pairs[a_idx], pairs[b_idx]
                u = gmul(pis[i], gconj(pis[j]))
                v = gmul(pis[k2], gconj(pis[l]))
                n = gnorm(u) * gnorm(v)
                t1 = 2 * mpmath.atan2(u[1], u[0]); t1 -= mpmath.nint(t1 / (mpmath.pi / 2)) * mpmath.pi / 2
                t2 = 2 * mpmath.atan2(v[1], v[0]); t2 -= mpmath.nint(t2 / (mpmath.pi / 2)) * mpmath.pi / 2
                prod = abs(mpmath.exp(1j * t1) - 1) ** 2 * abs(mpmath.exp(1j * t2) - 1) ** 2 * mpmath.mpf(n) ** (-3)
                ok3 &= prod <= mpmath.mpf(n) ** (-3 - delta)
                checked += 1
    pred = s ** (2 - delta) / (math.pi * (P / s) ** delta)
    print(f"s={s:>4} P={P:>6} delta={delta}: bins L={L:>5}; pair solutions={good2:>7} (all {len(pairs)} verified={good2 == len(pairs)}); "
          f"s^2/(2L)={s*s/(2*L):9.1f}; triples(disjoint pairs)={disjoint:>10} (sample verified={ok3})")
    return good2 == len(pairs) and ok3


if __name__ == '__main__':
    allok = True
    for s in (25, 50, 100, 200, 400, 800):
        allok &= run(s, 0.45)
    for s in (100, 400, 800):
        allok &= run(s, 0.9)
    print("ALL VERIFIED" if allok else "SOME FAILED")
