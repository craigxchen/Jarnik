"""Combinatorics for the S_M-isotypic decomposition of measures on shells / slices.

Conventions
-----------
* A partition lam of M is a weakly decreasing tuple.  Boxes are (r, c), 0-based.
* The fixed tableau t places the positions 0..M-1 in row-reading order:
  pos(r, c) = lam_0 + ... + lam_{r-1} + c.
* A standard Young tableau S is stored as a dict box -> entry (1..M).
* ATY index (Ariki-Terasoma-Yamada, "Higher Specht polynomials", 1997): read the columns of
  S from bottom to top, left to right, giving a word w.  index(1) = 0; index(k+1) = index(k)
  if k+1 stands to the right of k in w, else index(k) + 1.  i(S) = tableau of indices.
  charge(S) = sum of indices = degree of the higher Specht polynomial F_t^S.
* Shell S^M(p,n) = {k in Z^M : sum k^+ = p, sum k^- = n}; s = p - n, t = p + n.
  Slice Q^M(p,n) = disjoint union of shells S^M(p-i, n-i), 0 <= i <= min(p,n).
* An S_M-orbit is represented by its weakly decreasing value tuple v.
"""
from functools import lru_cache
from math import factorial
from collections import Counter
import itertools


def partitions(n, maxpart=None, maxlen=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    if maxlen == 0:
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k, None if maxlen is None else maxlen - 1):
            yield (k,) + rest


def conjugate(lam):
    if not lam:
        return ()
    return tuple(sum(1 for x in lam if x > j) for j in range(lam[0]))


def hook_dim(lam):
    n = sum(lam)
    conj = conjugate(lam)
    h = 1
    for i, r in enumerate(lam):
        for j in range(r):
            h *= (r - j - 1) + (conj[j] - i - 1) + 1
    return factorial(n) // h


def boxes(lam):
    return [(r, c) for r in range(len(lam)) for c in range(lam[r])]


def pos_of(lam):
    """dict box -> position under row reading."""
    out, k = {}, 0
    for r in range(len(lam)):
        for c in range(lam[r]):
            out[(r, c)] = k
            k += 1
    return out


def row_positions(lam):
    P = pos_of(lam)
    return [[P[(r, c)] for c in range(lam[r])] for r in range(len(lam))]


def col_positions(lam):
    P = pos_of(lam)
    conj = conjugate(lam)
    return [[P[(r, c)] for r in range(conj[c])] for c in range(len(conj))]


def standard_tableaux(lam):
    """all SYT of shape lam, as dicts box -> entry."""
    n = sum(lam)
    out = []

    def rec(shape, k, filling):
        if k > n:
            out.append(dict(filling))
            return
        for r in range(len(lam)):
            c = shape[r]
            if c < lam[r] and (r == 0 or shape[r - 1] > c):
                shape[r] += 1
                filling[(r, c)] = k
                rec(shape, k + 1, filling)
                del filling[(r, c)]
                shape[r] -= 1

    rec([0] * len(lam), 1, {})
    return out


def aty_index(lam, S):
    """index tableau i(S) (dict box -> index) and charge."""
    conj = conjugate(lam)
    word = []
    for c in range(len(conj)):
        for r in range(conj[c] - 1, -1, -1):
            word.append(S[(r, c)])
    where = {letter: i for i, letter in enumerate(word)}
    n = len(word)
    idx = {1: 0}
    for k in range(1, n):
        idx[k + 1] = idx[k] if where[k + 1] > where[k] else idx[k] + 1
    I = {b: idx[S[b]] for b in S}
    return I, sum(idx.values())


def major_index(lam, S):
    """maj(S) = sum of descents k (k+1 in a lower row than k)."""
    rowof = {S[b]: b[0] for b in S}
    n = len(S)
    return sum(k for k in range(1, n) if rowof[k + 1] > rowof[k])


def ssyt_fillings(lam, content):
    """semistandard tableaux of shape lam with given content.
    content: list of (value, multiplicity) with values strictly increasing.
    Returns a list of dicts box -> value (rows weakly increasing, columns strictly)."""
    L = len(lam)
    out = []

    def rec(i, shape, filling):
        if i == len(content):
            if list(shape) == list(lam):
                out.append(dict(filling))
            return
        val, m = content[i]
        # add a horizontal strip of size m to shape inside lam
        newshape = list(shape)

        def strip(r, rem):
            if r == L:
                if rem == 0:
                    added = []
                    for rr in range(L):
                        for c in range(shape[rr], newshape[rr]):
                            filling[(rr, c)] = val
                            added.append((rr, c))
                    rec(i + 1, tuple(newshape), filling)
                    for b in added:
                        del filling[b]
                return
            # row r can grow up to min(lam[r], shape[r-1]) (old shape of row above: horizontal strip)
            cap = lam[r] if r == 0 else min(lam[r], shape[r - 1])
            for a in range(0, min(rem, cap - shape[r]) + 1):
                newshape[r] = shape[r] + a
                strip(r + 1, rem - a)
            newshape[r] = shape[r]

        strip(0, m)

    rec(0, tuple([0] * L), {})
    return out


@lru_cache(None)
def kostka(lam, mu):
    """K_{lam,mu}, mu a tuple of multiplicities (any order)."""
    content = [(i, m) for i, m in enumerate(mu)]
    return len(ssyt_fillings(lam, content))


def shell_orbits(M, p, n):
    """weakly decreasing value tuples of the S_M-orbits of S^M(p,n)."""
    out = set()
    Ps = list(partitions(p, maxlen=M)) if p > 0 else [()]
    for P in Ps:
        Ns = list(partitions(n, maxlen=M - len(P))) if n > 0 else [()]
        for Nn in Ns:
            z = M - len(P) - len(Nn)
            if z < 0:
                continue
            v = tuple(sorted(list(P) + [-x for x in Nn] + [0] * z, reverse=True))
            out.add(v)
    return sorted(out, reverse=True)


def slice_orbits(M, p, n):
    out = []
    for i in range(0, min(p, n) + 1):
        out += shell_orbits(M, p - i, n - i)
    return out


def content_of(v):
    c = Counter(v)
    return [(val, c[val]) for val in sorted(c)]


def mu_of(v):
    return tuple(sorted(Counter(v).values(), reverse=True))


def orbit_size(v):
    c = Counter(v)
    s = factorial(len(v))
    for m in c.values():
        s //= factorial(m)
    return s


def shell_points(M, p, n):
    """all points of S^M(p,n) (for direct cross-checks)."""
    out = []
    for v in shell_orbits(M, p, n):
        out += sorted(set(itertools.permutations(v)))
    return out


def phi_floor(M, p, n):
    """floor of phi_M(p,n) = (M-1) min_{1<=j<=M-1} (p/(M-j) + n/j)  (upper.md, Theorem 3.1;
    used here only as a degree guess / comparison, never as an input to a claim)."""
    from fractions import Fraction
    if M == 1:
        return 0
    v = min(Fraction(p, M - j) + Fraction(n, j) for j in range(1, M))
    return int((M - 1) * v)
