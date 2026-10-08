"""Search integer cotangent cliques with a fixed scale L.

For each A in [L, Amax] factor S=A^2+L^2 by a quadratic sieve, list
divisors d <= dcap*A, candidates X=A+d, and enumerate all cliques of
size m (finite coordinates, including A) with least element A.  For each
clique compute the exact all-edge lcm N and report

  Cstar = 2 N^(1/4) arctan(L/A)         (normalised arc / sqrt(R))
  K5    = maxchord / N^(1/5)            (five-point CG scale, m=4)
  expo  = log N / log(A/L)

Usage: python3 search_fixedL.py L Amax m dcap
"""
import math
import sys
from itertools import combinations

import numpy as np

from cot_lib import all_edge_lcm, is_clique, edge_norm


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].tolist()


def sqrt_minus1_mod(p):
    # p = 1 mod 4: find r with r^2 = -1
    for a in range(2, p):
        r = pow(a, (p - 1) // 4, p)
        if r * r % p == p - 1:
            return r
    raise ValueError


def factor_table(L, Amax):
    """fac[A] = dict prime->exp for S = A^2+L^2, A in [0, Amax]."""
    S = [A * A + L * L for A in range(Amax + 1)]
    rem = S[:]
    fac = [dict() for _ in range(Amax + 1)]
    bound = math.isqrt(Amax * Amax + L * L) + 1
    for p in primes_upto(bound):
        if p == 2 or L % p == 0 or p % 4 == 1:
            # candidate residues A mod p with p | A^2+L^2
            if L % p == 0:
                roots = [0]
            elif p == 2:
                roots = [1 if L % 2 else 0]
            else:
                r = sqrt_minus1_mod(p)
                roots = sorted({(r * L) % p, (-r * L) % p})
            for r0 in roots:
                for A in range(r0, Amax + 1, p):
                    x = rem[A]
                    e = 0
                    while x % p == 0:
                        x //= p
                        e += 1
                    if e:
                        rem[A] = x
                        fac[A][p] = fac[A].get(p, 0) + e
    for A in range(Amax + 1):
        if rem[A] > 1:
            fac[A][rem[A]] = fac[A].get(rem[A], 0) + 1
    return fac


def divisors(f, cap):
    ds = [1]
    for p, e in f.items():
        new = []
        for d in ds:
            q = d
            for _ in range(e + 1):
                if q > cap:
                    break
                new.append(q)
                q *= p
        ds = new
    return ds


def maxchord_ratio(N, L, A):
    g = math.gcd(A, L)
    s = L // g
    n = edge_norm(A, L)
    a = A // g
    eps = 2 if (a % 2 and s % 2) else 1
    chord2 = 4 * N * s * s / (eps * n)
    return math.sqrt(chord2) / N ** 0.2


def search(L, Amax, m, dcap, report=15):
    fac = factor_table(L, Amax)
    res = []
    count = 0
    for A in range(max(L, 1), Amax + 1):
        ds = divisors(fac[A], dcap * A)
        cand = sorted(A + d for d in ds)
        # adjacency among candidates
        adj = {x: set() for x in cand}
        for x, y in combinations(cand, 2):
            if (x * y + L * L) % (y - x) == 0:
                adj[x].add(y)
                adj[y].add(x)
        # cliques of size m-1 among candidates
        def extend(cl, pool):
            if len(cl) == m - 1:
                yield cl
                return
            for y in sorted(pool):
                if cl and y <= cl[-1]:
                    continue
                yield from extend(cl + [y], pool & adj[y])
        for cl in extend([], set(cand)):
            xs = [A] + cl
            count += 1
            N = all_edge_lcm(xs, L)
            if N <= 1:
                continue
            Cs = 2 * N ** 0.25 * math.atan(L / A)
            K5 = maxchord_ratio(N, L, A)
            expo = math.log(N) / math.log(A / L) if A > L else float('inf')
            res.append((Cs, K5, expo, N, tuple(xs)))
    return res, count


if __name__ == '__main__':
    L, Amax, m, dcap = (int(a) for a in sys.argv[1:5])
    res, count = search(L, Amax, m, dcap)
    print('L=%d Amax=%d m=%d dcap=%d cliques=%d' % (L, Amax, m, dcap, count))
    print('-- smallest C_* --')
    for r in sorted(res)[:12]:
        print('C*=%.4f K5=%.4f expo=%.4f N=%d X=%s' % r)
    print('-- smallest K5 (with N>10^6) --')
    for r in sorted(res, key=lambda t: t[1]):
        if r[3] > 10 ** 6:
            print('C*=%.4f K5=%.4f expo=%.4f N=%d X=%s' % r)
            break
    big = [r for r in res if r[3] > 10 ** 8]
    for r in sorted(big, key=lambda t: t[1])[:12]:
        print('C*=%.4f K5=%.4f expo=%.4f N=%d X=%s' % r)
