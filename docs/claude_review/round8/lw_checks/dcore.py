"""Core library for lw.md checks: Paley skew-Hadamard matrices, the skew double core,
multi-basis profiles, and exact support computations (integer arithmetic only).

Index 0 is the bordered row/column ("infinity"); indices 1..q correspond to F_q elements 0..q-1.
"""
import itertools, random
from fractions import Fraction


def legendre(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1


def paley_skew(q):
    """Normalised skew-Hadamard K = I + S of order M = q+1 (q = 3 mod 4 prime), K(x,0) = 1."""
    assert q % 4 == 3
    M = q + 1
    S = [[0] * M for _ in range(M)]
    for x in range(1, M):
        S[x][0] = 1
        S[0][x] = -1
    for x in range(q):
        for y in range(q):
            S[x + 1][y + 1] = legendre(y - x, q)  # Jacobsthal Q(x,y) = chi(y-x), skew as chi(-1) = -1
    K = [[S[x][y] + (1 if x == y else 0) for y in range(M)] for x in range(M)]
    # checks: skew, Hadamard, first column ones
    for x in range(M):
        assert K[x][0] == 1
        for y in range(M):
            assert S[x][y] == -S[y][x]
    for x in range(M):
        for y in range(M):
            assert sum(K[x][k] * K[y][k] for k in range(M)) == (M if x == y else 0)
    return K


def double_core_types(K):
    """The 2(M-1) cut vectors (+-1) of exact.md Theorem 4.1: cuts of columns (a,0) of K and of K^T."""
    M = len(K)
    fam1 = [[K[x][a] * K[x][0] for x in range(M)] for a in range(1, M)]
    KT = [[K[y][x] for y in range(M)] for x in range(M)]
    fam2 = [[KT[x][a] * KT[x][0] for x in range(M)] for a in range(1, M)]
    return fam1, fam2


def image(c, types):
    """Integer image t_j = <c, s_j> for each type vector."""
    return [sum(ci * si for ci, si in zip(c, s)) for s in types]


def support_size(c, types):
    return sum(1 for t in image(c, types) if t != 0)


def rank_frac(rows):
    """Exact rank over Q of a list of integer rows (Fraction Gaussian elimination)."""
    A = [[Fraction(v) for v in r] for r in rows]
    rk = 0
    ncol = len(A[0]) if A else 0
    for col in range(ncol):
        piv = None
        for i in range(rk, len(A)):
            if A[i][col] != 0:
                piv = i
                break
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        for i in range(len(A)):
            if i != rk and A[i][col] != 0:
                f = A[i][col] / A[rk][col]
                A[i] = [a - f * b for a, b in zip(A[i], A[rk])]
        rk += 1
    return rk


def nullvec_int(rows, n):
    """A primitive integer vector spanning the kernel of `rows` (assumed 1-dimensional), or None."""
    A = [[Fraction(v) for v in r] for r in rows]
    pivcols = []
    rk = 0
    for col in range(n):
        piv = None
        for i in range(rk, len(A)):
            if A[i][col] != 0:
                piv = i
                break
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        p = A[rk][col]
        A[rk] = [a / p for a in A[rk]]
        for i in range(len(A)):
            if i != rk and A[i][col] != 0:
                f = A[i][col]
                A[i] = [a - f * b for a, b in zip(A[i], A[rk])]
        pivcols.append(col)
        rk += 1
    free = [c for c in range(n) if c not in pivcols]
    if len(free) != 1:
        return None
    fcol = free[0]
    v = [Fraction(0)] * n
    v[fcol] = Fraction(1)
    for i, pc in enumerate(pivcols):
        v[pc] = -A[i][fcol]
    from math import gcd
    den = 1
    for x in v:
        den = den * x.denominator // gcd(den, x.denominator)
    iv = [int(x * den) for x in v]
    g = 0
    for x in iv:
        g = gcd(g, abs(x))
    iv = [x // g for x in iv]
    # sign normalisation
    for x in iv:
        if x != 0:
            if x < 0:
                iv = [-y for y in iv]
            break
    return iv


def random_permuted_bases(K, t, seed):
    """t bases: basis 0 = nonconstant columns of K (as cuts with column 0); others random row permutations."""
    M = len(K)
    rng = random.Random(seed)
    bases = []
    for i in range(t):
        perm = list(range(M))
        if i > 0:
            rng.shuffle(perm)
        bases.append([[K[perm[x]][a] * K[perm[x]][0] for x in range(M)] for a in range(1, M)])
    return bases
