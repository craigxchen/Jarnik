"""Peeling recursion B_M(p,n) (upper bound for reg Q^M(p,n)) and the semigroup lower bound L_M(p,n).

B_1(p,n) = 0.
B_M(p,n) = max_i ( r_(i) + i - 1 ),  r_(1) >= r_(2) >= ... the sorted values of
           { B_{M-1}(p-c, n) : 0 <= c <= p }  U  { B_{M-1}(p, n-c) : 1 <= c <= n }.
(slices k_1 = c of Q^M(p,n); sorted-decreasing order is optimal for max(r_i + i - 1)).

Semigroup lower bound: generators (cost (T_b, T_a), ord C(a+b+1,2)), a+b+1 <= M (antisymmetrised
orbit of (-a..b) on a+b+1 coordinates), plus (1,0;1), (0,1;1).  L_M(p,n) = max total ord with
total cost <= (p,n) componentwise.
"""
from functools import lru_cache
from fractions import Fraction as F
import sys

@lru_cache(maxsize=None)
def B(M, p, n):
    if M == 1:
        return 0
    vals = [B(M - 1, p - c, n) for c in range(0, p + 1)] + [B(M - 1, p, n - c) for c in range(1, n + 1)]
    vals.sort(reverse=True)
    return max(v + i for i, v in enumerate(vals))

def T(x):
    return x * (x + 1) // 2

def gens(M):
    g = []
    for r in range(2, M + 1):
        for a in range(0, r):
            b = r - 1 - a
            g.append((T(b), T(a), r * (r - 1) // 2))
    return g

@lru_cache(maxsize=None)
def L(M, p, n):
    best = 0
    for (cp, cn, o) in gens(M):
        if cp <= p and cn <= n and (cp, cn) != (0, 0):
            best = max(best, o + L(M, p - cp, n - cn))
    return best

def phi(M, p, n):
    """LP bound: max sum lambda_g ord_g s.t. sum lambda_g cost_g <= (p,n); uses <= 2 generators."""
    g = gens(M)
    best = F(0)
    for (a1, b1, o1) in g:
        # single generator
        lam = min([F(p, a1)] if a1 else [] + [F(n, b1)] if b1 else [])  # placeholder
    # do it properly: enumerate pairs
    best = F(0)
    for i, (a1, b1, o1) in enumerate(g):
        lims = []
        if a1: lims.append(F(p, a1))
        if b1: lims.append(F(n, b1))
        lam = min(lims)
        best = max(best, lam * o1)
        for (a2, b2, o2) in g[i + 1:]:
            det = a1 * b2 - a2 * b1
            if det == 0:
                continue
            l1 = F(p * b2 - n * a2, det)
            l2 = F(a1 * n - b1 * p, det)
            if l1 >= 0 and l2 >= 0:
                best = max(best, l1 * o1 + l2 * o2)
    return best

if __name__ == '__main__':
    sys.setrecursionlimit(100000)
    M = int(sys.argv[1]); tmax = int(sys.argv[2])
    for t in range(1, tmax + 1):
        row = []
        for n in range(0, t // 2 + 1):
            p = t - n
            row.append(f"({p},{n}):B={B(M,p,n)},L={L(M,p,n)},phi={float(phi(M,p,n)):.2f}")
        print(f"t={t}: " + "  ".join(row))
