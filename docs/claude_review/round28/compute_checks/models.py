"""Profiles (pattern sets) for slab boxes, record conventions (round 11 bigger.py / round 27 killer k6_zono.py).
A profile is a sorted tuple of canonical subsets S of {0..n-1}; the box uses |kfull(S)| <= L for every S.
Restrictions to a point subset J (a subcluster's profile): S -> S cap J, re-indexed in the order of J (the last
element of J becomes the base point), trivial (empty/full) columns dropped, canonicalised, deduplicated."""
import itertools, random
import numpy as np
from slib import canon

def legendre(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q - 1) // 2, q) == 1 else -1)

def paley_K(q):
    idx = ['inf'] + list(range(q)); M = q + 1
    K = np.zeros((M, M), dtype=int)
    for a, x in enumerate(idx):
        for b, y in enumerate(idx):
            if x == 'inf': K[a, b] = 1
            elif y == 'inf': K[a, b] = -1
            else: K[a, b] = (1 if x == y else 0) + legendre(y - x, q)
    if not ((K + K.T) == 2 * np.eye(M, dtype=int)).all():
        K = K.T
    assert ((K + K.T) == 2 * np.eye(M, dtype=int)).all()
    return K

def sdc(q):
    """skew double core (round 12 fake.md; killer k6_zono.sdc7 for q = 7)"""
    K = paley_K(q); M = q + 1
    Dg = np.diag(K[:, 0]); K = Dg @ K @ Dg
    pats = set()
    for a in range(1, M):
        S = [x for x in range(M) if K[x, a] == -1]
        pats.add(canon(M, tuple(S))); pats.add(canon(M, tuple(sorted(set(S) ^ {0, a}))))
    return M, tuple(sorted(pats))

def had_paley(q):
    """Paley-I skew Hadamard matrix of order q+1 (= paley_K), normalised to first row 1; sign columns = rows 2.."""
    K = paley_K(q); M = q + 1
    H = K * K[0][None, :]          # multiply columns so that the first row is all ones
    assert (H[0] == 1).all() and ((H @ H.T) == M * np.eye(M, dtype=int)).all()
    pats = set(canon(M, tuple(x for x in range(M) if H[r, x] == 1)) for r in range(1, M))
    return M, tuple(sorted(pats))

def syl(k):
    H = [[1]]
    for _ in range(k):
        H = [r + r for r in H] + [r + [-x for x in r] for r in H]
    M = 2 ** k
    return M, tuple(sorted(set(canon(M, tuple(x for x in range(M) if r[x] == 1)) for r in H[1:])))

def restrict(M, pats, J):
    """profile of the subcluster J (list of point indices; order fixes the new indexing, last = base)"""
    pos = {x: i for i, x in enumerate(J)}
    n = len(J)
    out = set()
    for S in pats:
        T = tuple(sorted(pos[x] for x in S if x in pos))
        if 0 < len(T) < n:
            out.add(canon(n, T))
    return n, tuple(sorted(out))

def random_profile(n, seed):
    """round-27 killer law (k8_random.py): r ~ U{n-1..2n} distinct canonical patterns, each point in S with prob 1/2
    (nonempty proper), L ~ U{1,2}.  Returns (r, L, pats)."""
    rng = random.Random(seed)
    r = rng.randint(n - 1, 2 * n)
    pats = set()
    while len(pats) < r:
        S = tuple(x for x in range(n) if rng.random() < 0.5)
        if 0 < len(S) < n:
            pats.add(canon(n, S))
    L = rng.choice([1, 2])
    return r, L, tuple(sorted(pats))

def fullrank(n, pats):
    """rank of the profile vectors (indicator columns plus the all-ones vector) = n ?"""
    import sympy
    rows = [[1] * n] + [[1 if x in S else 0 for x in range(n)] for S in pats]
    return sympy.Matrix(rows).rank() == n
