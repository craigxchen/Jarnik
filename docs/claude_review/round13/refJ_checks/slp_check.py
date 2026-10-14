"""Referee checks for J6(b).
(1) SLP of B = Q[x_1..x_s]/(x_i^3) with l = sum eps_i x_i: rank of x l^t : B_i -> B_{i+t} equals min(h_i, h_{i+t})
    for all i, t (exact, s <= 6, several sign vectors).  (Proof in refJ.md is the sl_2 argument; this is a sanity check.)
(2) Hilbert function of R/J(A) for one-cut boxes A = {k in {-1,0,1}^d : |sum_{i in S} eps_i k_i| <= L} and for
    disjoint multi-cut boxes, computed exactly from evaluation ranks, equals the product of Lefschetz truncations
    max(h_i - h_{i-2L-1}, 0) (cut blocks) and (1+t+t^2) (free coordinates).
(3) strict unimodality of trinomial coefficients h^{(s)}_i for s >= 2, i < s (s <= 200), used for rho = s + L.
"""
import itertools, math
import numpy as np
from flint import fmpz_mat
out = []
def P(s): print(s); out.append(s)

def trinom(s):
    h = [1]
    for _ in range(s):
        h = [sum(h[i - t] for t in range(3) if 0 <= i - t < len(h)) for i in range(len(h) + 2)]
    return h

def monos(nv, i):
    res = []
    for c in itertools.combinations_with_replacement(range(nv), i):
        e = [0] * nv
        for x in c: e[x] += 1
        res.append(tuple(e))
    return res

# (1) SLP
def slp_check(s, eps):
    basis = {i: [e for e in monos(s, i) if max(e) <= 2] for i in range(2 * s + 1)}
    h = [len(basis[i]) for i in range(2 * s + 1)]
    assert h == trinom(s)
    ok = True
    for t in range(1, 2 * s + 1):
        for i in range(0, 2 * s + 1 - t):
            # matrix of l^t : B_i -> B_{i+t}
            idx = {e: r for r, e in enumerate(basis[i + t])}
            cols = []
            for e in basis[i]:
                poly = {e: 1}
                for _ in range(t):
                    new = {}
                    for mono, c in poly.items():
                        for x in range(s):
                            m2 = list(mono); m2[x] += 1
                            if m2[x] <= 2:
                                new[tuple(m2)] = new.get(tuple(m2), 0) + c * eps[x]
                    poly = {k: v for k, v in new.items() if v}
                col = [0] * len(basis[i + t])
                for mono, c in poly.items(): col[idx[mono]] = c
                cols.append(col)
            if not cols or not basis[i + t]: continue
            M = fmpz_mat(len(cols), len(basis[i + t]), [x for c in cols for x in c])
            if M.rank() != min(h[i], h[i + t]): ok = False
    return ok

for s in range(1, 7):
    for eps in [tuple([1] * s), tuple([1, -1] * s)[:s], tuple([-1] * s)]:
        P(f"SLP s={s} eps={eps}: max rank for all (i,t): {slp_check(s, eps)}")

# (2) Hilbert functions
def hilbert_from_eval(A, d, maxdeg):
    Aa = np.array(A, dtype=np.int64); cols = []; ranks = []
    for i in range(maxdeg + 1):
        for e in monos(d, i):
            if max(e) <= 2:
                col = np.ones(len(A), dtype=np.int64)
                for x, ex in enumerate(e):
                    if ex: col = col * Aa[:, x] ** ex
                cols.append(col)
        M = np.array(cols, dtype=np.int64)
        ranks.append(fmpz_mat(M.shape[0], M.shape[1], [int(t) for t in M.flatten()]).rank())
    return [ranks[0]] + [ranks[i] - ranks[i - 1] for i in range(1, maxdeg + 1)]

def polymul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i + j] += x * y
    return r

def lefschetz(s, L):
    h = trinom(s); N = 2 * L + 1
    return [max(h[i] - (h[i - N] if i >= N else 0), 0) for i in range(len(h))]

def predicted(blocks, free):
    p = [1]
    for (s, L) in blocks: p = polymul(p, lefschetz(s, L))
    for _ in range(free): p = polymul(p, [1, 1, 1])
    while len(p) > 1 and p[-1] == 0: p.pop()
    return p

cases = [  # (d, list of (S, eps, L))
    (5, [((0, 1, 2, 3, 4), (1, 1, 1, 1, 1), 1)]),
    (5, [((0, 1, 2, 3, 4), (1, -1, 1, -1, -1), 2)]),
    (6, [((0, 1, 2, 3), (1, -1, -1, 1), 1)]),
    (6, [((0, 1, 2), (1, 1, -1), 0), ((3, 4, 5), (1, -1, 1), 1)]),
    (7, [((0, 1, 2), (1, 1, 1), 1), ((3, 4, 5, 6), (-1, 1, -1, 1), 1)]),
    (7, [((0, 1), (1, -1), 0), ((2, 3, 4), (1, 1, 1), 1), ((5, 6), (1, 1), 1)]),
    (7, [((0, 1, 2, 3, 4, 5), (1, -1, 1, -1, 1, -1), 2)]),
]
for d, cuts in cases:
    A = [k for k in itertools.product((-1, 0, 1), repeat=d)
         if all(abs(sum(e * k[i] for i, e in zip(S, eps))) <= L for S, eps, L in cuts)]
    Hc = hilbert_from_eval(A, d, 2 * d)
    while len(Hc) > 1 and Hc[-1] == 0: Hc.pop()
    free = d - sum(len(S) for S, _, _ in cuts)
    Hp = predicted([(len(S), L) for S, _, L in cuts], free)
    P(f"d={d} cuts={[(S, eps, L) for S, eps, L in cuts]} |A|={len(A)}: H(R/J(A)) exact = {Hc}; Lefschetz product = {Hp}; equal: {Hc == Hp}")

# (3) strict unimodality
bad = []
for s in range(2, 201):
    h = trinom(s)
    if not all(h[i] < h[i + 1] for i in range(s)): bad.append(s)
P(f"strict increase h_i < h_(i+1) for 0 <= i < s, all 2 <= s <= 200: {not bad}")
open(__file__.replace('.py', '_out.txt'), 'w').write("\n".join(out) + "\n")
