"""Core exact library for the round-11 finite test of ZE.

Conventions (match round10/determinant.md):
  n points, d = n-1, k' = (k_1..k_d) in Z^d, k_n = -sum k'.
  Box B_T = {|k_x| <= T, x < n}.  A sign column beta in {+-1}^n gives <k,beta> = 2 k(S), S = {x: beta_x=+1}.
  A slab is (S, L): |k(S)| <= L  (equivalently |<k,beta>| <= w = 2L).
  h_A(beta) = 2 * (max_A k(S) - min_A k(S)).
  V(A,m) = {c in Z^A : sum_k c_k P(k') = 0 for all P of degree < m}.
"""
import itertools, math
from functools import lru_cache
import flint

def kfull(kp):
    return tuple(kp) + (-sum(kp),)

def kS(kp, S):
    kf = kfull(kp)
    return sum(kf[x] for x in S)

def build_A(n, T, slabs):
    """slabs: list of (S, L) with S a tuple of point indices in 0..n-1 (0-based)."""
    d = n - 1
    A = []
    for kp in itertools.product(range(-T, T + 1), repeat=d):
        if all(abs(kS(kp, S)) <= L for S, L in slabs):
            A.append(kp)
    return A

def monomials(d, deg):
    """exponent tuples alpha in N^d with |alpha| == deg"""
    out = []
    for c in itertools.combinations_with_replacement(range(d), deg):
        a = [0] * d
        for i in c:
            a[i] += 1
        out.append(tuple(a))
    return sorted(set(out), reverse=True)

def monomials_lt(d, m):
    out = []
    for t in range(m):
        out += monomials(d, t)
    return out

def mono_val(k, a):
    v = 1
    for ki, ai in zip(k, a):
        if ai:
            v *= ki ** ai
    return v

def moment_matrix(A, m):
    d = len(A[0])
    rows = monomials_lt(d, m)
    return flint.fmpz_mat([[mono_val(k, a) for k in A] for a in rows])

def hilbert_A(A, mmax):
    """H_A(j) = dim Poly_{<= j}|_A for j = 0..mmax (exact rank)."""
    out = []
    for j in range(mmax + 1):
        out.append(moment_matrix(A, j + 1).rank())
    return out

def kernel_basis(A, m, reduce=True):
    """exact Z-basis of the integer kernel lattice V(A,m) cap Z^A (saturated), LLL-reduced.
    Method: unimodular row reduction (HNF) of [M^T | I_U]; rows with zero left block span ker_Z M."""
    M = moment_matrix(A, m)
    E, U = M.nrows(), M.ncols()
    rows = []
    for j in range(U):
        rows.append([int(M[i, j]) for i in range(E)] + [1 if t == j else 0 for t in range(U)])
    H = flint.fmpz_mat(rows).hnf()
    ker = []
    for i in range(H.nrows()):
        if all(H[i, t] == 0 for t in range(E)):
            v = [int(H[i, E + t]) for t in range(U)]
            if any(v):
                ker.append(v)
    if not ker:
        return []
    B = flint.fmpz_mat(ker)
    if reduce:
        B = B.lll()
    out = [[int(B[i, j]) for j in range(U)] for i in range(B.nrows())]
    assert len(out) == U - M.rank()
    return out

def check_kernel(A, m, basis):
    M = moment_matrix(A, m)
    for c in basis:
        v = M * flint.fmpz_mat([[x] for x in c])
        assert all(v[i, 0] == 0 for i in range(v.nrows()))
    return True

def width(A, S):
    vals = [kS(k, S) for k in A]
    return max(vals) - min(vals)

def siegel_H(A, T, m):
    U = len(A)
    return (U * T ** (m - 1)) ** 2

def phi_forms(A, m, basis, d):
    """Coefficient vectors (in monomials of degree m) of Phi_c(u) = sum_k c_k <k,u>^m for c in basis.
    Returns (mons, list of integer coefficient lists) spanning the tangent-cone linear system (row-reduced)."""
    mons = monomials(d, m)
    multin = []
    for a in mons:
        v = math.factorial(m)
        for ai in a:
            v //= math.factorial(ai)
        multin.append(v)
    Mm = flint.fmpz_mat([[mono_val(k, a) for k in A] for a in mons])
    C = flint.fmpz_mat([[x] for x in basis[0]]) if basis else None
    cols = []
    for c in basis:
        v = Mm * flint.fmpz_mat([[x] for x in c])
        cols.append([multin[i] * int(v[i, 0]) for i in range(len(mons))])
    if not cols:
        return mons, []
    F = flint.fmpz_mat(cols)
    r = F.rank()
    # integer row basis via HNF-like LLL then drop zero rows
    R = F.lll()
    rows = [[int(R[i, j]) for j in range(R.ncols())] for i in range(R.nrows())]
    rows = [rw for rw in rows if any(rw)]
    assert len(rows) == r, (len(rows), r)
    return mons, rows
