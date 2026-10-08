"""Exact check of Lemma 5.2 (anchoring) in ptolemy.md.

For random configurations of primitive integer vectors w_0=(1,0), w_1..w_m (arbitrary, with
deliberately shared p-adic structure), put D = lcm_j |y_j| (the choice r_p = max_j v_p(y_j)) and
P_i = primitive part of (D x_i, y_i).  Check, at every prime p dividing some determinant:
  (a) P_i represent the same point of M_{0,k} (all cross-ratios unchanged);
  (b) det(P_0,P_j) = +-1;
  (c) v_p(det(P_i,P_j)) = delta_p + sum_{T contains i,j} l_T(p)  with one constant delta_p >= 0,
      where l_T(p) are the tree edge lengths (boundary splits, side T not containing 0).
"""
import random, itertools, sys
from fractions import Fraction
from math import gcd

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 5)

def vp(n, p):
    n = abs(n); v = 0
    while n % p == 0:
        n //= p; v += 1
    return v

def primes_of(n):
    n = abs(n); out = set(); q = 2
    while q * q <= n:
        while n % q == 0:
            out.add(q); n //= q
        q += 1
    if n > 1:
        out.add(n)
    return out

def lcm(a, b):
    return a * b // gcd(a, b)

def det(u, v):
    return u[0] * v[1] - u[1] * v[0]

def edge_lengths(g, k):
    rows = list(range(1, k))
    out = {}
    for r in range(2, k - 1):
        for T in itertools.combinations(rows, r):
            Tc = [x for x in range(k) if x not in T]
            best = None
            for a, b in itertools.combinations(T, 2):
                for c, d in itertools.combinations(Tc, 2):
                    I = min(g[a, b] + g[c, d] - g[a, c] - g[b, d], g[a, b] + g[c, d] - g[a, d] - g[b, c])
                    best = I if best is None else min(best, I)
            out[frozenset(T)] = max(0, best)
    return out

checks = 0
for trial in range(300):
    k = random.randint(5, 7)
    base = [2, 3, 5, 7, 11, 13]
    w = [(1, 0)]
    while len(w) < k:
        # build vectors with shared prime-power structure: x = sum of products of small primes
        x = random.choice([1, -1]) * random.randint(1, 10**4) * random.choice([1, 2, 3, 6, 5, 25, 7])
        y = random.randint(1, 10**3) * random.choice([1, 4, 9, 27, 5, 49])
        g0 = gcd(x, y)
        v = (x // g0, y // g0)
        if all(det(v, u) != 0 for u in w):
            w.append(v)
    D = 1
    for v in w[1:]:
        D = lcm(D, abs(v[1]))
    P = [(1, 0)]
    for (x, y) in w[1:]:
        a, b = D * x, y
        g0 = gcd(a, b)
        P.append((a // g0, b // g0))
    # (a) cross-ratios unchanged
    for a, b, c, d in itertools.combinations(range(k), 4):
        cr = lambda V: Fraction(det(V[a], V[b]) * det(V[c], V[d]), det(V[a], V[c]) * det(V[b], V[d]))
        assert cr(w) == cr(P)
    # (b)
    for j in range(1, k):
        assert abs(det(P[0], P[j])) == 1
    # (c)
    primes = set()
    for i, j in itertools.combinations(range(k), 2):
        primes |= primes_of(det(P[i], P[j]))
    for p in primes:
        g = {}
        for i, j in itertools.combinations(range(k), 2):
            g[i, j] = g[j, i] = vp(det(P[i], P[j]), p)
        ell = edge_lengths(g, k)
        deltas = set()
        for i, j in itertools.combinations(range(1, k), 2):
            s = sum(l for T, l in ell.items() if i in T and j in T)
            deltas.add(g[i, j] - s)
        assert len(deltas) == 1 and min(deltas) >= 0, (p, deltas)
        checks += 1
print(f"Lemma 5.2 verified: 300 configurations, {checks} prime checks; det(P_0,P_j) = +-1 always")
