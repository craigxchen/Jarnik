"""Exact rational-function reconstruction of a univariate function given by exact evaluations."""
from fractions import Fraction as Fr
import qpoly as Q

def nullspace(M):
    """exact nullspace of matrix M (list of rows of Fractions)."""
    M = [row[:] for row in M]
    nr, nc = len(M), len(M[0])
    piv = []
    r = 0
    for c in range(nc):
        p = None
        for i in range(r, nr):
            if M[i][c] != 0:
                p = i; break
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        inv = 1 / M[r][c]
        M[r] = [x * inv for x in M[r]]
        for i in range(nr):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == nr: break
    free = [c for c in range(nc) if c not in piv]
    basis = []
    for f in free:
        v = [Fr(0)] * nc
        v[f] = Fr(1)
        for i, c in enumerate(piv):
            v[c] = -M[i][f]
        basis.append(v)
    return basis

def reconstruct(f, maxdeg=40, start=3):
    """find p/q (coprime, q monic) with f(x) = p(x)/q(x). Tries total degree D = m + n increasing."""
    xs = [Fr(k, 7) + 3 for k in range(2 * maxdeg + 10)]
    vals = {}
    def F(x):
        if x not in vals: vals[x] = Fr(f(x))
        return vals[x]
    for D in range(0, maxdeg + 1):
        for n in range(0, D + 1):
            m = D - n
            pts = xs[:m + n + 4]
            rows = []
            for x in pts:
                y = F(x)
                rows.append([x ** i for i in range(m + 1)] + [-y * x ** j for j in range(n + 1)])
            ns = nullspace(rows)
            if len(ns) == 1:
                v = ns[0]
                p = Q.norm(v[:m + 1]); q = Q.norm(v[m + 1:])
                if q:
                    g = Q.gcd(p, q) if p else q
                    lc = q[-1]
                    p = Q.scal(1 / lc, p); q = Q.scal(1 / lc, q)
                    return p, q
    return None
