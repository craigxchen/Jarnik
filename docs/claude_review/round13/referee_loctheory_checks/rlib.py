"""Referee library (independent of the author's lfi.py).
Exact leading-form spaces J(X)_e of finite point sets X in Z^b, pullbacks, ideal degree pieces modulo a
monomial (grid) ideal, ranks exact (fmpz) or mod p (lower bounds on rational rank)."""
import itertools
from math import comb
import flint

P = 2147483647  # prime 2^31-1


def monos_eq(b, e, cap=None):
    """exponent tuples of total degree e in b vars, each exponent < cap[q] if cap given"""
    out = []
    def rec(q, left, cur):
        if q == b - 1:
            if cap is None or left < cap[q]:
                out.append(tuple(cur + [left]))
            return
        hi = left if cap is None else min(left, cap[q] - 1)
        for a in range(hi, -1, -1):
            rec(q + 1, left - a, cur + [a])
    if b == 0:
        return [()] if e == 0 else []
    rec(0, e, [])
    return out


def monos_le(b, e, cap=None):
    out = []
    for d in range(e + 1):
        out += monos_eq(b, d, cap)
    return out


def evalmono(pt, ex):
    v = 1
    for x, a in zip(pt, ex):
        if a:
            v *= x ** a
    return v


def h_func_exact(X, e, cap=None):
    """h_X(e) = rank of evaluation matrix of monomials of degree <= e on X (exact over Q)."""
    ms = monos_le(len(X[0]), e, cap)
    M = flint.fmpz_mat([[evalmono(p, m) for m in ms] for p in X])
    return M.rank()


def h_func_modp(X, e, cap=None, p=P):
    ms = monos_le(len(X[0]), e, cap)
    M = flint.nmod_mat([[evalmono(pt, m) % p for m in ms] for pt in X], p)
    return M.rank()


def J_space(X, e, cap=None):
    """basis (list of dicts exponent->int) of J(X)_e modulo nothing: top-degree parts of polys of degree <= e
    vanishing on X.  If cap given, only reduced monomials are used (valid when X lies in the grid with
    sizes cap, since then the grid polynomials reduce every function); result is J(X)_e mod grid."""
    b = len(X[0])
    ms = monos_le(b, e, cap)
    M = flint.fmpz_mat([[evalmono(p, m) for m in ms] for p in X])
    N, nul = M.nullspace()
    top = [i for i, m in enumerate(ms) if sum(m) == e]
    if nul == 0 or not top:
        return []
    # rows of nullspace vectors restricted to top monomials
    rows = [[N[i, c] for i in top] for c in range(nul)]
    R = flint.fmpz_mat(rows)
    Rr, den, rk = R.rref()
    basis = []
    for r in range(rk):
        d = {}
        den = None
        vals = [Rr[r, c] for c in range(len(top))]
        # rref over fmpz returns (matrix, den) style? handle fmpq
        for c, v in enumerate(vals):
            if v != 0:
                d[ms[top[c]]] = v
        basis.append(d)
    return basis


def poly_mul(f, g):
    out = {}
    for a, x in f.items():
        for b_, y in g.items():
            k = tuple(i + j for i, j in zip(a, b_))
            out[k] = out.get(k, 0) + x * y
    return {k: v for k, v in out.items() if v != 0}


def linform_pow_cache():
    return {}


def pullback(f, dirs, b):
    """f: dict over exponents in len(dirs) vars; dirs: list of integer vectors of length b.
    returns dict in b vars (homogeneous if f is)."""
    lin = []
    for v in dirs:
        lin.append({tuple(1 if i == q else 0 for i in range(b)): int(v[q]) for q in range(b) if v[q] != 0})
    out = {}
    one = {tuple([0] * b): 1}
    powcache = {}
    def lp(q, a):
        if (q, a) not in powcache:
            if a == 0:
                powcache[(q, a)] = one
            else:
                powcache[(q, a)] = poly_mul(lp(q, a - 1), lin[q])
        return powcache[(q, a)]
    for ex, c in f.items():
        term = {tuple([0] * b): c}
        for q, a in enumerate(ex):
            if a:
                term = poly_mul(term, lp(q, a))
        for k, v in term.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}


def reduce_cap(f, cap):
    return {k: v for k, v in f.items() if all(a < c for a, c in zip(k, cap))}


def ideal_degree_rows(gens_by_deg, b, D, cap):
    """gens_by_deg: dict e -> list of homogeneous dicts (in b vars).  Returns list of rows (dicts) spanning
    (ideal generated)_D modulo the monomial ideal of exponents >= cap."""
    rows = []
    for e, gl in gens_by_deg.items():
        if e > D:
            continue
        mults = monos_eq(b, D - e, cap)
        for g in gl:
            for m in mults:
                r = reduce_cap({tuple(i + j for i, j in zip(k, m)): v for k, v in g.items()}, cap)
                if r:
                    rows.append(r)
    return rows


def rank_rows(rows, cols, modp=True, p=P):
    idx = {c: i for i, c in enumerate(cols)}
    if not rows:
        return 0
    if modp:
        M = flint.nmod_mat(len(rows), len(cols), p)
        for i, r in enumerate(rows):
            for k, v in r.items():
                M[i, idx[k]] = v % p
        return M.rank()
    M = flint.fmpz_mat(len(rows), len(cols))
    for i, r in enumerate(rows):
        for k, v in r.items():
            M[i, idx[k]] = v
    return M.rank()
