"""Exact small T = 1 models: Hilbert function, V(A,m) basis, tangent forms.

model(name, L) -> dict with A, d, H (Hilbert function of S/J, exact over Q via fmpq rref when |A| small,
else mod p), and helper to get an exact integer basis of V(A,m).
"""
import sys, json, time
import numpy as np
import flint
from flib import *

P = 2305843009213693951


def profile(name):
    if name.startswith('sdc'):
        return sdc_patterns(int(name[3:]))
    if name.startswith('pal'):
        return pal_patterns(int(name[3:]))
    if name.startswith('syl'):
        return sylvester_patterns({4: 2, 8: 3, 16: 4}[int(name[3:])])
    raise ValueError(name)


def graded_monomials(d, maxdeg):
    mons = []
    for j in range(maxdeg + 1):
        mons += reduced_monomials(d, j)
    return mons


def evalmat(A, mons):
    Aa = np.array(A, dtype=np.int64)
    cols = []
    for a in mons:
        v = np.ones(len(A), dtype=np.int64)
        for x, e in enumerate(a):
            if e:
                v = v * Aa[:, x] ** e
        cols.append(v)
    return np.stack(cols, axis=1)


def hilbert_modp(A, d, p=P, exact=False):
    """Hilbert function via rank profile of the graded evaluation matrix (columns in degree order)."""
    j = 0
    H = []
    prev = 0
    while prev < len(A):
        mons = graded_monomials(d, j)
        E = evalmat(A, mons)
        if exact:
            M = flint.fmpz_mat(E.tolist())
            r = M.rank()
        else:
            M = flint.nmod_mat(E.shape[0], E.shape[1], [int(x) % p for x in E.ravel()], p)
            r = M.rank()
        H.append(r - prev)
        prev = r
        j += 1
    return H


def V_basis(A, d, m):
    """exact integer basis (LLL-reduced) of V(A,m) = {c : sum_k c_k k^a = 0 for all reduced a, |a| < m}."""
    mons = graded_monomials(d, m - 1)
    E = evalmat(A, mons)
    Mt = flint.fmpz_mat(E.T.tolist())  # E x |A|
    X, nul = Mt.nullspace()
    U = len(A)
    if nul == 0:
        return []
    B = flint.fmpz_mat([[X[i, jj] for i in range(U)] for jj in range(nul)])
    # make primitive integer rows, then LLL
    B = B.lll()
    basis = [[int(B[i, jj]) for jj in range(U)] for i in range(B.nrows())]
    Ev = np.array(basis, dtype=object).dot(np.array(E, dtype=object))
    assert all(x == 0 for x in np.array(Ev).ravel())
    return basis


if __name__ == '__main__':
    name = sys.argv[1]; L = int(sys.argv[2])
    M, pats = profile(name)
    A = box(M, pats, L)
    d = M - 1
    t = time.time()
    H = hilbert_modp(A, d, exact=(len(A) <= 2500))
    print(name, 'L=%d |A|=%d H=%s rho=%d (%.1fs)' % (L, len(A), H, len(H) - 1, time.time() - t))
