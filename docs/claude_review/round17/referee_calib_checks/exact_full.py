#!/usr/bin/env python3
"""Referee: fully exact, symmetry-free computation of Psi_A on the Fourier circle C_1.

No eigenspace decomposition, no multimodular step.  For a box A (list of k in {-1,0,1}^q):
 * M = monomial evaluation matrix (k^alpha, alpha in {0,1,2}^q, |alpha| < m), integer;
 * exact rank and pivot columns S from rref of the integer Gram M^T M (python-flint, exact);
 * W_{i,j}(k) = coefficient of zeta^j in F_i(k) = khat_1^i khat_{-1}^{m-i}, computed in the group ring Z[Z_q];
 * residual Gram Gam = W^T W - (M_S^T W)^T (M_S^T M_S)^{-1} (M_S^T W), exact rationals (fmpq);
 * <R_i, R_i'> = sum_{j,j'} Gam[(i,j),(i',j')] zeta^{j-j'} reduced in Q(zeta) (= Q[Z_q]/(1+z+...+z^{q-1})).
Checks: <R_i,R_i'> = 0 in Q(zeta) for i != i' (exact theta-independence); <R_i,R_i> rational;
Psi_A(C_1) = 4^{-m} sum binom(m,i)^2 ||R_i||^2 / U  (probability normalisation; ||.||^2 here is the plain sum).
"""
import sys, os, time, itertools
from math import comb
from fractions import Fraction
import numpy as np
import flint

def monomials(q, m):
    out = []
    def rec(x, left, cur):
        if x == q:
            out.append(tuple(cur)); return
        for e in range(0, min(2, left) + 1):
            cur.append(e); rec(x + 1, left - e, cur); cur.pop()
    for deg in range(m):
        tmp = []
        out_before = len(out)
        rec(0, deg, [])
    # rec above collects |alpha| <= deg; dedupe
    s = sorted(set(a for a in out if sum(a) < m), key=lambda a: (sum(a), a))
    return s

def mono_matrix(A, mons):
    A = A.astype(np.int64)
    M = np.ones((A.shape[0], len(mons)), dtype=np.int64)
    for c, al in enumerate(mons):
        col = np.ones(A.shape[0], dtype=np.int64)
        for x, e in enumerate(al):
            if e == 1: col *= A[:, x]
            elif e == 2: col *= A[:, x] ** 2
        M[:, c] = col
    return M

def groupring_powers(A, q, m):
    """W[i] : array (U, q) integer, coefficients of F_i(k) = khat_1^i khat_{-1}^{m-i} in Z[Z_q]."""
    U = A.shape[0]
    A = A.astype(object)
    k1 = np.zeros((U, q), dtype=object); km1 = np.zeros((U, q), dtype=object)
    for x in range(q):
        k1[:, x] += A[:, x]; km1[:, (-x) % q] += A[:, x]
    def mul(a, b):
        c = np.zeros_like(a)
        for s in range(q):
            c += a[:, [s]] * np.roll(b, s, axis=1)  # (a_s z^s) * b
        return c
    one = np.zeros((U, q), dtype=object); one[:, 0] = 1
    p1 = [one]; pm = [one]
    for e in range(m):
        p1.append(mul(p1[-1], k1)); pm.append(mul(pm[-1], km1))
    return [mul(p1[i], pm[m - i]) for i in range(m + 1)]

def reduce_cyclo(c, q):
    """c: list of q Fractions = sum c_t z^t in Q[Z_q]; return (is_rational, value) in Q(zeta_q)."""
    rest = c[1:]
    if all(r == rest[0] for r in rest):
        return True, c[0] - rest[0]
    return False, None

def run(A, q, m, log):
    t0 = time.time()
    U = A.shape[0]
    mons = monomials(q, m)
    M = mono_matrix(A, mons)
    Gfull = flint.fmpz_mat((M.T @ M).tolist())
    R, den, _rk = Gfull.rref()
    rank = 0; piv = []
    nr, nc = R.nrows(), R.ncols()
    row = 0
    for c in range(nc):
        if row < nr and R[row, c] != 0:
            piv.append(c); row += 1
    rank = len(piv)
    MS = M[:, piv]
    W = groupring_powers(A, q, m)
    Wcat = np.hstack([w for w in W]).astype(object)  # U x ((m+1) q)
    # exact integer products
    MSo = MS.astype(object)
    G = flint.fmpz_mat((MSo.T.dot(MSo)).tolist()) if rank else None
    B = MSo.T.dot(Wcat) if rank else None
    WW = Wcat.T.dot(Wcat)
    nW = Wcat.shape[1]
    if rank:
        X = flint.fmpq_mat(G).solve(flint.fmpq_mat(flint.fmpz_mat(B.tolist())))
        # verify G X == B exactly
        assert flint.fmpq_mat(G) * X == flint.fmpq_mat(flint.fmpz_mat(B.tolist()))
        BtX = flint.fmpq_mat(flint.fmpz_mat(B.T.tolist())) * X
    def gam(a, b):
        v = Fraction(int(WW[a, b]))
        if rank:
            e = BtX[a, b]
            v -= Fraction(int(e.p), int(e.q))
        return v
    def ip(i, ip_, resid=True):
        c = [Fraction(0)] * q
        for j in range(q):
            for jp in range(q):
                a, b = i * q + j, ip_ * q + jp
                v = gam(a, b) if resid else Fraction(int(WW[a, b]))
                c[(j - jp) % q] += v
        return reduce_cyclo(c, q)
    out = {'U': U, 'rank': rank, 'dimV': U - rank, 'nmon': len(mons)}
    offdiag_ok = True
    for i in range(m + 1):
        for i2 in range(i + 1, m + 1):
            isr, val = ip(i, i2)
            if not (isr and val == 0):
                offdiag_ok = False
    Rn = []; Fn = []
    for i in range(m + 1):
        isr, val = ip(i, i); assert isr, 'diag not rational'
        Rn.append(val)
        isr, val = ip(i, i, resid=False); assert isr
        Fn.append(val)
    Psi = sum(comb(m, i) ** 2 * Rn[i] for i in range(m + 1)) / (Fraction(4) ** m * U)
    E2m = sum(comb(m, i) ** 2 * Fn[i] for i in range(m + 1)) / (Fraction(4) ** m * U)
    out.update(offdiag_zero=offdiag_ok, Psi=Psi, E2m=E2m, Rn=[r / U for r in Rn])
    log(f'q={q} m={m}: U={U} nmon={len(mons)} rank={rank} dimV={U-rank} offdiag<R_i,R_i\'>=0 exactly: {offdiag_ok}')
    log(f'   ||R_i||^2 (prob. norm) = ' + ', '.join(f'{r}' for r in out['Rn']))
    log(f'   Psi_A(C_1) = {Psi} ~ {float(Psi):.10e};  E l^2m = {float(E2m):.10e}; rel = {float(Psi/E2m):.6e}  ({time.time()-t0:.1f}s)')
    return out

if __name__ == '__main__':
    q, L = int(sys.argv[1]), int(sys.argv[2])
    ms = [int(x) for x in sys.argv[3].split(',')]
    here = os.path.dirname(os.path.abspath(__file__))
    if L == 1 and not os.path.exists(os.path.join(here, 'data', f'A_q{q}_L1.npy')):
        rows = [[0] * q]
        for x in range(q):
            for s in (1, -1):
                e = [0] * q; e[x] = s; rows.append(e)
            for y in range(q):
                if y != x:
                    e = [0] * q; e[x] = 1; e[y] = -1; rows.append(e)
        A = np.array(rows, dtype=np.int8)
    else:
        A = np.load(os.path.join(here, 'data', f'A_q{q}_L{L}.npy'))
    f = open(os.path.join(here, 'out', f'exact_full_q{q}_L{L}.txt'), 'a')
    def log(s):
        print(s, flush=True); f.write(s + '\n'); f.flush()
    for m in ms:
        run(A, q, m, log)
