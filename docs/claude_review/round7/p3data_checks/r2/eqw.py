"""Equal-weight code search for the prime-floor problem W_sc (round 7, p3data, finishing session).

For a weight tau >= w0 and r columns, an equal-weight profile S in {+-1}^(M x r) has margins
mu_xy = tau (d_xy - r/2), so mu >= D  <=>  d_xy >= T_xy(r) := ceil(r/2 + D_xy/tau).
* r_P(tau): least r with sum_{x<y} T_xy(r) <= r floor(M/2) ceil(M/2) (Plotkin: a column separates at
  most floor(M/2) ceil(M/2) pairs).  Exact lower bound for equal-weight profiles of weight tau.
* SA (sacode.c) searches S for r = r_P, r_P + 1, ...; the first success is an upper bound
  r tau for W_sc (equal weights tau >= w0 are admissible).
* LP polish: for the found S, the best weights >= w0 (fls.lp_S, slack must be 0); value <= r tau.
Every S is re-checked in Python with exact integer distances against T (float only in D/tau).
"""
import os, sys, math, ctypes
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from fls import lp_S                      # noqa: E402

_so = os.path.join(HERE, 'libsacode.so')
if not os.path.exists(_so) or os.path.getmtime(_so) < os.path.getmtime(os.path.join(HERE, 'sacode.c')):
    os.system('gcc -O3 -march=native -shared -fPIC -o %s %s -lm' % (_so, os.path.join(HERE, 'sacode.c')))
_lib = ctypes.CDLL(_so)
_lib.sa.restype = ctypes.c_long
_lib.sa.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int8),
                    ctypes.c_long, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_uint64]


def thresholds(D, tau, r):
    T = np.ceil(r / 2.0 + D / tau - 1e-12).astype(np.int32)
    np.fill_diagonal(T, 0)
    return T


def plotkin_r(D, tau, rmax=100000):
    M = D.shape[0]
    iu = np.triu_indices(M, 1)
    cap = (M // 2) * ((M + 1) // 2)
    for r in range(1, rmax):
        if thresholds(D, tau, r)[iu].astype(np.int64).sum() <= r * cap:
            return r
    return None


def random_balanced(M, r, rng):
    S = np.empty((M, r), dtype=np.int8)
    half = np.array([1] * (M // 2) + [-1] * (M - M // 2), dtype=np.int8)
    for j in range(r):
        S[:, j] = rng.permutation(half)
    return S


def sa_run(T, S, iters, T0=1.0, T1=0.03, pf=0.5, seed=1):
    M, r = S.shape
    Tc = np.ascontiguousarray(T, dtype=np.int32)
    Sc = np.ascontiguousarray(S, dtype=np.int8)
    v = _lib.sa(M, r, Tc.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                Sc.ctypes.data_as(ctypes.POINTER(ctypes.c_int8)), int(iters), T0, T1, pf, seed)
    return int(v), Sc


def exact_ok(S, T):
    M = S.shape[0]
    d = (S[:, None, :] != S[None, :, :]).sum(axis=2)
    iu = np.triu_indices(M, 1)
    cols_ok = all(abs(int(S[:, j].astype(int).sum())) < M for j in range(S.shape[1]))
    return bool((d[iu] >= T[iu]).all()) and cols_ok


def search_tau(D, tau, rng, iters=3_000_000, restarts=3, extra=40, seed0=0):
    """least r (from r_P up) for which SA finds a code; returns (r_P, r, S)."""
    M = D.shape[0]
    rP = plotkin_r(D, tau)
    S_prev = None
    for r in range(rP, rP + extra + 1):
        T = thresholds(D, tau, r)
        for k in range(restarts):
            if S_prev is not None and k == 0:      # warm start: previous best plus one random column
                S0 = np.hstack([S_prev, random_balanced(M, 1, rng)])
            else:
                S0 = random_balanced(M, r, rng)
            v, S = sa_run(T, S0, iters, seed=seed0 * 7919 + r * 101 + k)
            if v == 0 and exact_ok(S, T):
                return rP, r, S
            if k == 0:
                S_best = S
        S_prev = S_best
    return rP, None, None


def polish(D, S, w0):
    val, w, sl, _ = lp_S(D, S.astype(np.int64), w0)
    if w is None or sl > 1e-7:
        return None, None
    return float(val), w
