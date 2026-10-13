# Shared exact/numeric helpers for round-11 subcluster checks.
import itertools, math, numpy as np
from fractions import Fraction

def monomials(d, m):
    """exponent tuples alpha in N^d with |alpha| < m"""
    return [a for t in range(m) for a in itertools.product(range(t+1), repeat=d) if sum(a) == t]

def rank_mod_p(M, p=2147483629):
    """rank modulo a prime < 2^31 (vectorised numpy int64); rank mod p <= rank over Q."""
    A = np.array([[int(x) % p for x in row] for row in M], dtype=np.int64)
    nr, nc = A.shape; r = 0
    for c in range(nc):
        nzr = np.nonzero(A[r:, c])[0]
        if len(nzr) == 0: continue
        piv = r + nzr[0]
        if piv != r: A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p-2, p)
        A[r] = (A[r]*inv) % p
        f = A[:, c].copy(); f[r] = 0
        nz = np.nonzero(f)[0]
        if len(nz): A[nz] = (A[nz] - (f[nz, None]*A[r][None, :]) % p) % p
        r += 1
        if r == nr: break
    return r

def slab_box(J, base, cols, T, w):
    """J: list of point indices (len n); base in J; cols: list of 0/1 tuples over all points (global profile).
    Returns A as list of full k-vectors in Z^J_0 (ordered as J), with |k_x|<=T for x!=base and
    |<k,beta_p|J>| <= w for each column, beta=2a-1 (on the zero-sum slice <k,beta> = 2<k,a>)."""
    n = len(J); bi = J.index(base)
    others = [i for i in range(n) if i != bi]
    A = []
    for ks in itertools.product(range(-T, T+1), repeat=n-1):
        k = [0]*n
        for i, v in zip(others, ks): k[i] = v
        k[bi] = -sum(ks)
        ok = True
        for a in cols:
            s = 2*sum(k[i]*a[J[i]] for i in range(n))
            if abs(s) > w: ok = False; break
        if ok: A.append(tuple(k))
    return A

def width(A, vec):
    vals = [sum(ki*vi for ki, vi in zip(k, vec)) for k in A]
    return max(vals) - min(vals)

def moment_matrix(A, base_idx, m):
    """rows: monomials of degree < m in the coordinates k' (drop base coordinate); cols: points of A."""
    n = len(A[0]); others = [i for i in range(n) if i != base_idx]
    mons = monomials(n-1, m)
    Mat = [[math.prod(k[o]**e for o, e in zip(others, al)) for k in A] for al in mons]
    return Mat, mons

P31 = 2147483629

def rref_mod_p(A, p=P31):
    """A: int64 numpy array. Returns (R, pivots) reduced row echelon form mod p."""
    A = np.array(A, dtype=np.int64) % p
    nr, nc = A.shape; r = 0; piv = []
    for c in range(nc):
        if r == nr: break
        nzr = np.nonzero(A[r:, c])[0]
        if len(nzr) == 0: continue
        pr_ = r + nzr[0]
        if pr_ != r: A[[r, pr_]] = A[[pr_, r]]
        inv = pow(int(A[r, c]), p-2, p)
        A[r] = (A[r]*inv) % p
        f = A[:, c].copy(); f[r] = 0
        nz = np.nonzero(f)[0]
        if len(nz): A[nz] = (A[nz] - (f[nz, None]*A[r][None, :]) % p) % p
        piv.append(c); r += 1
    return A[:r], piv

def nullspace_mod_p(A, p=P31):
    """basis (rows) of {x : A x = 0} mod p."""
    A = np.array(A, dtype=np.int64)
    nc = A.shape[1]
    R, piv = rref_mod_p(A, p)
    free = [c for c in range(nc) if c not in set(piv)]
    N = np.zeros((len(free), nc), dtype=np.int64)
    for i, f in enumerate(free):
        N[i, f] = 1
        for j, pc in enumerate(piv):
            N[i, pc] = (-R[j, f]) % p
    return N

def rank_rows_mod_p(A, p=P31):
    if len(A) == 0: return 0
    R, piv = rref_mod_p(np.array(A, dtype=np.int64), p)
    return len(piv)
