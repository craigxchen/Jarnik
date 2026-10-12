#!/usr/bin/env python3
"""R2. Prop 2.2 (gcd lemma) and Lemma 2.1 (sub-cluster identity) on ARBITRARY subsets, exact integers.

stability.md checks Prop 2.2 only on subsets of short windows.  Prop 2.2 is claimed for every subset B
of every configuration (with B's own chord constant C_B).  Here B ranges over
  * random subsets of all lattice points of the circle (huge C_B),
  * random subsets of angular windows (small C_B),
  * all subsets of the 8 angularly closest points to random centres,
on circles whose norms carry content (2, primes 3 mod 4, repeated split primes).
Exact form of Prop 2.2 (multiply (2.2) by s^2, write log C_B = (1/2) log X - W/4, X = max |z_x - z_y|^2):
      Norm(G_B)^(s^2) * exp(8 F*(s)) * N^(s(s-2))  <=  X^(2 s (s-1)).
Exact form of Lemma 2.1, with u = z/G_B, cuts and levels computed on u (not on z):
      prod |z_x-z_y|^8 = N^(s(s-2)) Norm(G_B)^(s^2) prod_{cuts of u} p^((2S-s)^2) prod |beta_xy|^8.
Also: equality in (2.2) for s = 2 iff |beta|^2 = 2 (Remark 2.3(1)).
Usage: python3 r2_gcd_subsets.py [samples_per_circle] [seed]
"""
import sys
import random
from math import atan2, pi
from rlib import (factorize, lattice_points, gsub, gnorm, ggcd, ggcd_list, gdivexact, split_prime,
                  vpi, gconj, F8int)

NORMS = [
    5 * 13 * 17 * 29 * 37,             # squarefree
    5 ** 3 * 13 ** 2 * 17,             # prime powers
    2 * 5 ** 2 * 13 * 17 * 29,          # ramified content
    9 * 5 * 13 * 17 * 29,               # inert content
    4 * 49 * 5 ** 2 * 13 * 17,          # mixed content
    5 ** 4 * 13 ** 2 * 29,
    65 * 85 * 145 * 41,                 # 5^3 13 17 29 41
    5 * 13 * 17 * 29 * 37 * 41 * 53,   # many points
    2 * 9 * 121 * 5 * 13 * 17 * 37,
    5 ** 2 * 13 ** 2 * 17 ** 2 * 29,
]


def cuts_of(us, Nb):
    """list of (p, S) over cuts (p, tau) of the configuration us (gcd 1, norm Nb)."""
    out = []
    for p, e in factorize(Nb).items():
        if p % 4 != 1:
            continue
        pi_ = split_prime(p)
        lev = [vpi(u, pi_) for u in us]
        for tau in range(1, e + 1):
            out.append((p, sum(1 for a in lev if a < tau)))
    return out


def test_subset(Z, N, stats):
    s = len(Z)
    G = ggcd_list(Z)
    nG = gnorm(G)
    X = max(gnorm(gsub(Z[a], Z[b])) for a in range(s) for b in range(a))
    lhs = nG ** (s * s) * F8int(s) * N ** (s * (s - 2))
    rhs = X ** (2 * s * (s - 1))
    stats['gcd'] += 1
    if lhs > rhs:
        stats['gcd_fail'] += 1
        print("  GCD LEMMA FAILS", N, Z)
    if s == 2:
        g = ggcd(Z[0], Z[1])
        beta2 = gnorm(gsub(Z[0], Z[1])) // gnorm(g)
        if (lhs == rhs) != (beta2 == 2):
            stats['eq_fail'] += 1
        if lhs == rhs:
            stats['equality_cases'] += 1
    # identity
    U = [gdivexact(z, G) for z in Z]
    Nb = N // nG
    L = 1
    R = N ** (s * (s - 2)) * nG ** (s * s)
    for a in range(s):
        for b in range(a):
            d = gsub(Z[a], Z[b])
            L *= gnorm(d) ** 4
            g = ggcd(Z[a], Z[b])
            R *= (gnorm(d) // gnorm(g)) ** 4
    for p, S in cuts_of(U, Nb):
        R *= p ** ((2 * S - s) ** 2)
    stats['ident'] += 1
    if L != R:
        stats['ident_fail'] += 1
        print("  IDENTITY FAILS", N, Z)


def main():
    samples = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    rng = random.Random(seed)
    stats = dict(gcd=0, gcd_fail=0, ident=0, ident_fail=0, eq_fail=0, equality_cases=0)
    for N in NORMS:
        fac = factorize(N)
        pts = sorted(lattice_points(fac), key=lambda z: atan2(z[1], z[0]) % (2 * pi))
        n = len(pts)
        # random subsets of the whole circle
        for _ in range(samples):
            s = rng.randint(2, min(9, n))
            test_subset(rng.sample(pts, s), N, stats)
        # random subsets of windows
        for _ in range(samples):
            w = rng.randint(3, min(14, n))
            i = rng.randrange(n)
            win = [pts[(i + t) % n] for t in range(w)]
            s = rng.randint(2, w)
            test_subset(rng.sample(win, s), N, stats)
        # all subsets of the 8 nearest points around a few centres
        from itertools import combinations
        for _ in range(3):
            i = rng.randrange(n)
            win = [pts[(i + t) % n] for t in range(min(8, n))]
            for s in range(2, len(win) + 1):
                for B in combinations(win, s):
                    test_subset(list(B), N, stats)
        print(f"N = {N} = {fac}: {n} points; cumulative {stats}")
    ok = stats['gcd_fail'] == 0 and stats['ident_fail'] == 0 and stats['eq_fail'] == 0
    print("ALL CHECKS PASSED:", ok)


if __name__ == "__main__":
    main()
