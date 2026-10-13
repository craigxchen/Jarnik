"""Core helpers for the discrete problem (Q) (round 7 p3data).

Cuts: nonconstant s in {+-1}^M modulo s ~ -s, encoded by a bitmask with bit 0 = 0 (vertex 0 on side +1);
bit x = 1 means s_x = -1.  Pair coefficient of cut s at pair (x,y): [s_x != s_y] - 1/2.
"""
import ctypes, os, math
import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
_lib = ctypes.CDLL(os.path.join(_here, 'libmaxcut.so'))
_lib.best_cuts.restype = ctypes.c_int
_lib.best_cuts.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.c_int,
                           ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ctypes.c_double)]


def best_cuts(nu, K=1):
    """Exact top-K weighted max-cuts of symmetric nonnegative nu (M <= 30 in practice)."""
    M = nu.shape[0]
    a = np.ascontiguousarray(nu, dtype=np.float64)
    masks = (ctypes.c_uint64 * K)()
    vals = (ctypes.c_double * K)()
    n = _lib.best_cuts(M, a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), K, masks, vals)
    return [int(masks[i]) for i in range(n)], [float(vals[i]) for i in range(n)]


def mask_to_s(mask, M):
    return np.array([-1 if (mask >> x) & 1 else 1 for x in range(M)], dtype=np.int64)


def s_to_mask(s):
    s = np.asarray(s)
    if s[0] < 0:
        s = -s
    return int(sum(1 << x for x in range(len(s)) if s[x] < 0))


def all_cut_masks(M):
    return list(range(2, 2 ** M, 2))   # bit 0 = 0, nonzero


def pair_index(M):
    return np.triu_indices(M, 1)


def pair_coeffs(masks, M):
    """Matrix P[pair, cut] = [s_x != s_y] - 1/2 (pairs in triu order)."""
    iu = pair_index(M)
    S = np.array([mask_to_s(m, M) for m in masks], dtype=np.int64)    # cuts x M
    return np.where(S[:, iu[0]] != S[:, iu[1]], 0.5, -0.5).T, S


def margins(masks, lam, M):
    """mu_xy = sum_s lam_s ([s_x != s_y] - 1/2) as an M x M matrix."""
    P, _ = pair_coeffs(masks, M)
    mu = np.zeros((M, M))
    iu = pair_index(M)
    v = P @ np.asarray(lam, dtype=float)
    mu[iu] = v
    mu[(iu[1], iu[0])] = v
    return mu


def primes_1mod4_above(X, count):
    out, n = [], X + 1
    while len(out) < count:
        if n % 4 == 1 and all(n % d for d in range(2, int(n ** 0.5) + 1)):
            out.append(n)
        n += 1
    return out
