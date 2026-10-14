"""T = 1 slab-box models.  A model = (d, cols, Ls): A = {k in {-1,0,1}^d : |<col_p, k>| <= L_p}.
Coordinates are the non-base points; the base coordinate k_base = -sum k is unconstrained except through columns."""
import itertools, math
import numpy as np


def cube(d):
    return np.array(list(itertools.product((-1, 0, 1), repeat=d)), dtype=np.int8)


def box(d, cols, Ls, chunk=None):
    C = np.array(cols, dtype=np.int64)
    Lv = np.array(Ls, dtype=np.int64)
    if 3 ** d <= 2_000_000:
        K = cube(d).astype(np.int64)
        ok = (np.abs(K @ C.T) <= Lv).all(axis=1)
        return K[ok]
    # chunked over the first 4 coordinates
    out = []
    rest = cube(d - 4).astype(np.int64)
    for head in itertools.product((-1, 0, 1), repeat=4):
        K = np.hstack([np.tile(np.array(head, dtype=np.int64), (len(rest), 1)), rest])
        ok = (np.abs(K @ C.T) <= Lv).all(axis=1)
        out.append(K[ok])
    return np.vstack(out)


def sylvester(k, L):
    """Syl(2^k, L): points F_2^k, base 0; coordinates x != 0 (as integers 1..n-1);
    columns k(S_a), S_a = {x : a.x = 1} (does not contain the base), |k(S_a)| <= L."""
    n = 2 ** k
    pts = list(range(1, n))
    cols = []
    for a in range(1, n):
        cols.append([1 if bin(a & x).count('1') % 2 == 1 else 0 for x in pts])
    return n - 1, cols, [L] * (n - 1)


def paley12(L):
    """Paley-12 (zetc §5.5): coordinates Z_11 (base = infinity), columns i + D, D = {0} u QR_11, |k(i+D)| <= L."""
    q = 11; D = [0, 1, 3, 4, 5, 9]
    cols = []
    for i in range(q):
        c = [0] * q
        for x in D:
            c[(i + x) % q] = 1
        cols.append(c)
    return q, cols, [L] * q


def cuts(d, sets, Ls):
    cols = [[1 if x in S else 0 for x in range(d)] for S in sets]
    return d, cols, list(Ls)
