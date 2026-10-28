"""Core exact routines for the round-18 census.

Conventions (ptolemy.md §1, §7; L8.md §1):
  k = m+1 points, anchor P_0 = (1,0) (= infinity), P_i = (X_i, Y_i) primitive integer vectors.
  Visible blocks: T subset of [m] = {1..m}, 2 <= |T| <= m-1  (2^m - m - 1 of them).
  Contact integer d_T(Phi) = prod_p p^{l_T(p)}, l_T(p) = p-adic tree edge length (Thm B(iv)).

In the integral chart (all Y_i = 1, which Lemma 5.2 always allows) the tree at p, rooted
towards the anchor, has an edge for T exactly at the depths t >= 1 for which T is a full residue
class of {X_i mod p^t}.  Hence
  l_T(p) = #{ t >= 1 : T is one residue class of (X_i mod p^t) },  2 <= |T| <= m-1,
and v_p(X_i - X_j) = delta_p + sum_{T contains i,j} l_T(p), delta_p = v_p(gcd of all differences)
(the invisible top block).  This is the exact Boolean form with t_ij = +-1, Y_j = 1.

Everything here is exact integer arithmetic; logs are used only for reporting.
"""
from itertools import combinations
from math import gcd, log
from flint import fmpz


def visible_blocks(m):
    """All visible T (as frozensets of labels 1..m), 2 <= |T| <= m-1."""
    out = []
    for s in range(2, m):
        for T in combinations(range(1, m + 1), s):
            out.append(frozenset(T))
    return out


def factor_primes(n):
    """Set of primes dividing the nonzero integer n (exact, via flint)."""
    n = abs(int(n))
    if n <= 1:
        return []
    return [int(p) for p, e in fmpz(n).factor()]


def vp(n, p):
    n = abs(int(n))
    if n == 0:
        raise ValueError("v_p(0)")
    c = 0
    while n % p == 0:
        n //= p
        c += 1
    return c


def contacts_integral(X, primes=None):
    """Contact profile of the configuration (inf, X_1..X_m), X distinct integers.

    Returns (l, delta) where l[p][T] = l_T(p) > 0 (visible T only) and delta[p] = top part.
    X is a list of length m (labels 1..m).  primes: optional list of primes to use
    (default: all primes dividing some difference).
    """
    m = len(X)
    if len(set(X)) != m:
        raise ValueError("points not distinct")
    if primes is None:
        ps = set()
        for i in range(m):
            for j in range(i + 1, m):
                ps.update(factor_primes(X[i] - X[j]))
        primes = sorted(ps)
    l = {}
    delta = {}
    for p in primes:
        vmax = max(vp(X[i] - X[j], p) for i in range(m) for j in range(i + 1, m))
        lp = {}
        dp = 0
        q = 1
        for t in range(1, vmax + 1):
            q *= p
            classes = {}
            for i, x in enumerate(X):
                classes.setdefault(x % q, []).append(i + 1)
            for cl in classes.values():
                if len(cl) == m:
                    dp += 1
                elif len(cl) >= 2:
                    T = frozenset(cl)
                    lp[T] = lp.get(T, 0) + 1
        if lp:
            l[p] = lp
        if dp:
            delta[p] = dp
    return l, delta


def contact_integers(X, m=None):
    """d_T for all visible T (dict), and the top part delta (int)."""
    m = len(X)
    l, delta = contacts_integral(X)
    d = {T: 1 for T in visible_blocks(m)}
    for p, lp in l.items():
        for T, e in lp.items():
            d[T] *= p ** e
    top = 1
    for p, e in delta.items():
        top *= p ** e
    return d, top, l


def check_boolean_identity(X):
    """Verify exactly X_i - X_j = +- top * prod_{T contains i,j} d_T (Lemma 5.2 / round4 L7 Thm 1)."""
    m = len(X)
    d, top, _ = contact_integers(X)
    for i in range(1, m + 1):
        for j in range(i + 1, m + 1):
            G = top
            for T, v in d.items():
                if i in T and j in T:
                    G *= v
            if abs(X[i - 1] - X[j - 1]) != G:
                return False
    return True


def admissible_at_prime(Ts, m):
    """Two-chain condition at one prime: the visible blocks Ts (with the top block [m] allowed)
    lie in at most two nested chains whose supports (unions) are disjoint."""
    Ts = sorted(set(Ts), key=len)
    # laminar is automatic for a p-adic tree; check chain decomposition
    # maximal elements of the laminar family:
    maxi = [T for T in Ts if not any(T < U for U in Ts)]
    if len(maxi) > 2:
        return False
    for M in maxi:
        inside = [T for T in Ts if T <= M]
        inside.sort(key=len)
        for a, b in zip(inside, inside[1:]):
            if not a < b:
                return False
    return True


def balance_stats(d):
    """min, max of log d_T over visible T (None if some d_T = 1)."""
    vals = [log(v) for v in d.values()]
    return min(vals), max(vals)


def kheight_k8(ell):
    """7*h_K finite part for k=8: sum over splits of (s(8-s)-14) * ell_T, s = |T| or |T^c| size.
    ell: dict T -> log d_T for visible T in [7]; split sizes: |T| and 8-|T|."""
    tot = 0.0
    for T, v in ell.items():
        s = min(len(T), 8 - len(T))
        tot += (s * (8 - s) - 14) * v
    return tot


def contact_integers_gcd(X):
    """Factoring-free contact integers in the integral chart (anchor infinity):
    g_T := gcd_{i,j in T}(X_i - X_j) = delta * prod_{T' visible, T' ⊇ T} d_{T'}, inverted top-down.
    Returns (d, delta)."""
    m = len(X)
    delta = 0
    for i in range(m):
        for j in range(i + 1, m):
            delta = gcd(delta, X[i] - X[j])
    delta = abs(delta)
    Ts = sorted(visible_blocks(m), key=len, reverse=True)
    d = {}
    for T in Ts:
        g = 0
        idx = sorted(T)
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                g = gcd(g, X[idx[a] - 1] - X[idx[b] - 1])
        sup = delta
        for T2, v in d.items():
            if T2 > T:
                sup *= v
        assert g % sup == 0
        d[T] = abs(g) // sup
    return d, delta


def integral_chart(points, anchor_index=0):
    """points: list of distinct Fractions (None = infinity).  Moves points[anchor_index] to infinity,
    returns the primitive integer representation X of the remaining points (gcd of differences 1)."""
    from fractions import Fraction as Fr
    a = points[anchor_index]
    us = []
    for i, z in enumerate(points):
        if i == anchor_index:
            continue
        if a is None:
            us.append(Fr(z))
        elif z is None:
            us.append(Fr(0))
        else:
            us.append(1 / (Fr(z) - Fr(a)))
    L = 1
    for u in us:
        L = L * u.denominator // gcd(L, u.denominator)
    X = [int(u * L) for u in us]
    g = 0
    for i in range(len(X)):
        for j in range(i + 1, len(X)):
            g = gcd(g, X[i] - X[j])
    return [x // g for x in X]
