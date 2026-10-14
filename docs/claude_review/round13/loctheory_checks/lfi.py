"""Leading-form-ideal utilities (exact over Q via python-flint fmpz_mat; optional mod-p).

For a finite X in Z^D:
  h_X(i)   = dim of {polys of degree <= i} restricted to X  (= rank of evaluation matrix)
  H_X(i)   = h_X(i) - h_X(i-1) = dim (S/J(X))_i
  J(X)_i   = top-degree parts of kernel vectors of the degree-<=i evaluation matrix.
All ranks are exact (fmpz) unless mod p is requested.
"""
import itertools
from functools import lru_cache
import flint

def monos(D, deg, cap=None):
    """exponent tuples of total degree exactly deg (each exponent <= cap if cap given)."""
    out = []
    def rec(pos, left, cur):
        if pos == D - 1:
            if cap is None or left <= cap:
                out.append(tuple(cur + [left]))
            return
        hi = left if cap is None else min(left, cap)
        for e in range(hi, -1, -1):
            rec(pos + 1, left - e, cur + [e])
    if D == 0:
        return [()] if deg == 0 else []
    rec(0, deg, [])
    return out

def evalrow(pt, mons):
    r = []
    for e in mons:
        v = 1
        for a, k in zip(pt, e):
            if k:
                v *= a ** k
        r.append(v)
    return r

def rank_exact(rows):
    if not rows or not rows[0]:
        return 0
    return flint.fmpz_mat(rows).rank()

def rank_modp(rows, p=2147475331):
    if not rows or not rows[0]:
        return 0
    return flint.nmod_mat(rows, p).rank()

def hilbert(X, maxdeg, cap=None, modp=False):
    """H_X(0..maxdeg). cap: exponent cap valid when X lies in a set where y^(cap+1) is a leading form
    (e.g. ternary cube cap=2); then monomials with larger exponents are not needed."""
    X = list(X)
    D = len(X[0])
    mons = []
    h = []
    rk = rank_modp if modp else rank_exact
    for i in range(maxdeg + 1):
        mons += monos(D, i, cap)
        rows = [evalrow(pt, mons) for pt in X]
        # rank of |X| x |mons|; transpose if wide is irrelevant for flint
        h.append(rk(rows))
        if h[-1] == len(X):
            break
    H = [h[0]] + [h[i] - h[i - 1] for i in range(1, len(h))]
    return H

def rho(X, cap=None, modp=False):
    H = hilbert(X, 10 ** 6, cap, modp)
    return len(H) - 1

def leading_forms(X, i, cap=None):
    """Basis (list of integer vectors over monos(D,i,cap)) of J(X)_i, exact."""
    X = list(X)
    D = len(X[0])
    allm = []
    for t in range(i + 1):
        allm += monos(D, t, cap)
    top = monos(D, i, cap)
    ntop = len(top)
    M = flint.fmpz_mat([evalrow(pt, allm) for pt in X])
    # kernel of M (columns = monomials): use fmpq nullspace via rref
    K, nul = M.nullspace()
    vecs = []
    nlow = len(allm) - ntop
    for c in range(nul):
        v = [int(K[r, c]) for r in range(len(allm))]
        vecs.append(v[nlow:])
    if not vecs:
        return top, []
    # reduce to a basis
    Mt = flint.fmpz_mat(vecs)
    r = Mt.rank()
    # take rref rows
    R, rr = flint.fmpq_mat(Mt).rref()
    basis = []
    for a in range(rr):
        row = [R[a, b] for b in range(ntop)]
        den = 1
        for q in row:
            den = den * int(q.q) // __import__('math').gcd(den, int(q.q))
        basis.append([int(q * den) for q in row])
    return top, basis

def ideal_dim_in_degree(i, D, gens_by_degree, cap=None):
    """dim of degree-i part of the ideal generated (in k[y], or in k[y]/(y^(cap+1)) if cap) by
    homogeneous generators given as dicts {exponent tuple: coeff}.  Exact."""
    top = monos(D, i, cap)
    idx = {e: a for a, e in enumerate(top)}
    rows = []
    for g in gens_by_degree:
        dg = sum(next(iter(g)))
        if dg > i:
            continue
        for mu in monos(D, i - dg, None):
            row = [0] * len(top)
            ok = False
            for e, c in g.items():
                f = tuple(a + b for a, b in zip(e, mu))
                if cap is not None and max(f) > cap:
                    continue
                row[idx[f]] += c
                ok = True
            if ok and any(row):
                rows.append(row)
    if not rows:
        return 0
    return rank_exact(rows)

P_DEFAULT = 2147475331

def leading_forms_modp(X, i, cap=None, p=P_DEFAULT):
    """mod-p basis of (leading forms of mod-p kernel) in degree i.  NOTE: the mod-p kernel can be LARGER than the
    rational one; use only inside a sandwich argument (see write-up)."""
    X = list(X)
    D = len(X[0])
    allm = []
    for t in range(i + 1):
        allm += monos(D, t, cap)
    top = monos(D, i, cap)
    ntop = len(top)
    M = flint.nmod_mat([[v % p for v in evalrow(pt, allm)] for pt in X], p)
    K, nul = M.nullspace()
    nlow = len(allm) - ntop
    vecs = [[int(K[r, c]) for r in range(nlow, len(allm))] for c in range(nul)]
    vecs = [v for v in vecs if any(v)]
    if not vecs:
        return top, []
    R, rr = flint.nmod_mat(vecs, p).rref()
    return top, [[int(R[a, b]) for b in range(ntop)] for a in range(rr)]

def ideal_dim_modp(i, D, gens, cap=None, p=P_DEFAULT):
    top = monos(D, i, cap)
    idx = {e: a for a, e in enumerate(top)}
    rows = []
    for g in gens:
        dg = sum(next(iter(g)))
        if dg > i:
            continue
        for mu in monos(D, i - dg, None):
            row = [0] * len(top)
            for e, c in g.items():
                f = tuple(a + b for a, b in zip(e, mu))
                if cap is not None and max(f) > cap:
                    continue
                row[idx[f]] = (row[idx[f]] + c) % p
            if any(row):
                rows.append(row)
    if not rows:
        return 0
    return flint.nmod_mat(rows, p).rank()
