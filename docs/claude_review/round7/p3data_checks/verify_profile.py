"""Exact verification of a candidate solution of the discrete problem (Q) (round 7, p3data).

A solution is a list of columns (cut mask, set of distinct primes): the column's sign pattern is the cut,
and each prime p in its set is one physical column of weight log p (e = 1), so W = sum log p.
Demands: integer unit-1 collision moduli m_xy (D_xy = 2 log m_xy).  C = sqrt 2, so c_C = 0.

Checks (all integer arithmetic except where marked):
 (P)  pair inequality  mu_xy >= 2 log m_xy  <=>  prod_sep p >= m^4 prod_same p.
 (S)  delta placement: rational t in [0,1]^M with, for every pair,
        a^4 (prod_sep p - m^2) >= m^4 prod_same p den^4,   |t_x - t_y| = a/den.
      This implies |sin((delta_x - delta_y)/2)| >= 2^(-1/2) m e^(-d/2) for delta = Delta t,
      Delta = C e^(-W/4) (proof in p3data.md section 3).
 (C4) every ternary zero-sum c of support <= 4 and every c = 2e_x - e_y - e_z has
      s_c = w(v(c)) - W/2 >= 0  <=>  prod p^(2|v_j|) >= N   (exact).
 (Call) min over all nonzero zero-sum integer c with |c_x| <= 30 of s_c, by MILP (floating point;
      reported as evidence-level certification).
Also: distinct primes = 1 mod 4 above X, nonconstant columns, distinct rows.
"""
import itertools, math
import numpy as np
from fractions import Fraction
from scipy.optimize import milp, LinearConstraint, Bounds
from gauss7 import is_prime
from qcore import mask_to_s


def expand(cols, M):
    """cols: list of (mask, [primes]) -> sign matrix S (M x r) and prime list"""
    S, ps = [], []
    for mask, plist in cols:
        s = mask_to_s(mask, M)
        for p in plist:
            S.append(s); ps.append(int(p))
    return np.array(S, dtype=np.int64).T, ps


def basic_checks(cols, M, X):
    S, ps = expand(cols, M)
    assert len(set(ps)) == len(ps), 'primes not distinct'
    for p in ps:
        assert p % 4 == 1 and p > X and is_prime(p), p
    for j in range(S.shape[1]):
        assert len(set(S[:, j])) == 2, 'constant column'
    assert len({tuple(r) for r in S}) == M, 'rows not distinct'
    return S, ps


def prods(S, ps, x, y):
    sep, same = 1, 1
    for j, p in enumerate(ps):
        if S[x, j] != S[y, j]:
            sep *= p
        else:
            same *= p
    return sep, same


def check_pairs(cols, mod, X):
    M = mod.shape[0]
    S, ps = basic_checks(cols, M, X)
    worst = None
    for x in range(M):
        for y in range(x + 1, M):
            sep, same = prods(S, ps, x, y)
            m = int(mod[x, y])
            ok = sep >= m ** 4 * same
            if not ok:
                return False, (x, y)
            sl = 0.5 * (math.log(sep) - math.log(same)) - 2 * math.log(m)   # mu - D (float, report only)
            worst = sl if worst is None else min(worst, sl)
    return True, worst


def check_placement(cols, mod, X, t):
    """t: list of Fractions in [0,1]."""
    M = mod.shape[0]
    S, ps = basic_checks(cols, M, X)
    assert all(0 <= v <= 1 for v in t)
    for x in range(M):
        for y in range(x + 1, M):
            sep, same = prods(S, ps, x, y)
            m = int(mod[x, y])
            dt = abs(t[x] - t[y])
            a, den = dt.numerator, dt.denominator
            if not (a ** 4 * (sep - m * m) >= m ** 4 * same * den ** 4):
                return False, (x, y)
    return True, None


def char_check_small(cols, M, X):
    S, ps = basic_checks(cols, M, X)
    A = (1 - S) // 2               # a_xj in {0,1}
    N = 1
    for p in ps:
        N *= p
    cs = []
    for sup in itertools.combinations(range(M), 4):
        for pos in itertools.combinations(sup, 2):
            neg = [x for x in sup if x not in pos]
            if pos[0] != sup[0]:
                continue           # c ~ -c
            c = np.zeros(M, dtype=np.int64); c[list(pos)] = 1; c[neg] = -1
            cs.append(c)
    for x in range(M):
        for y, z in itertools.combinations([u for u in range(M) if u != x], 2):
            c = np.zeros(M, dtype=np.int64); c[x] = 2; c[y] = -1; c[z] = -1
            cs.append(c)
    C = np.array(cs)
    V = np.abs(C @ A)              # |v_j|
    worst = None
    for i in range(len(cs)):
        lhs = 1
        for j, p in enumerate(ps):
            if V[i, j]:
                lhs *= p ** int(2 * V[i, j])
        if lhs < N:
            return False, cs[i]
        sc = float(sum(V[i, j] * math.log(p) for j, p in enumerate(ps)) - 0.5 * sum(math.log(p) for p in ps))
        worst = sc if worst is None else min(worst, sc)
    return True, (len(cs), worst)


def char_min_milp(cols, M, bound=30, time_limit=120):
    """min_{c != 0, sum c = 0, |c_x| <= bound} s_c via M MILPs (c_x0 >= 1).  Floating point."""
    S, ps = expand(cols, M)
    A = (1 - S) // 2
    w = np.array([math.log(p) for p in ps])
    r = len(ps)
    W = w.sum()
    best, argbest = None, None
    for x0 in range(M):
        # vars: c (M ints), t (r cont) ; minimize w.t ; t_j >= +-(A^T c)_j ; sum c = 0 ; c_x0 >= 1
        nvar = M + r
        cost = np.concatenate([np.zeros(M), w])
        rows, lo, hi = [], [], []
        At = A.T.astype(float)
        rows.append(np.hstack([-At, np.eye(r)])); lo.append(np.zeros(r)); hi.append(np.full(r, np.inf))
        rows.append(np.hstack([At, np.eye(r)])); lo.append(np.zeros(r)); hi.append(np.full(r, np.inf))
        rows.append(np.concatenate([np.ones(M), np.zeros(r)])[None, :]); lo.append([0]); hi.append([0])
        lb = np.concatenate([np.full(M, -bound), np.zeros(r)]); ub = np.concatenate([np.full(M, bound), np.full(r, np.inf)])
        lb[x0] = 1
        res = milp(cost, constraints=[LinearConstraint(np.vstack(rows), np.concatenate(lo), np.concatenate(hi))],
                   integrality=np.concatenate([np.ones(M), np.zeros(r)]), bounds=Bounds(lb, ub),
                   options=dict(time_limit=time_limit, disp=False))
        if res.x is None:
            continue
        val = res.fun - W / 2
        if best is None or val < best:
            best, argbest = val, np.round(res.x[:M]).astype(int)
    return best, argbest


def placement_milp(gam, time_limit=120):
    """min span L s.t. exists t in [0, L] with |t_x - t_y| >= gam_xy (order binaries).  Returns (L, t)."""
    M = gam.shape[0]
    pairs = [(x, y) for x in range(M) for y in range(x + 1, M)]
    nP = len(pairs)
    # vars: t (M), L, o (nP)
    nv = M + 1 + nP
    cost = np.zeros(nv); cost[M] = 1
    rows, lo, hi = [], [], []
    G = 2 + 2 * gam.max()
    for k, (x, y) in enumerate(pairs):
        g = gam[x, y]
        r1 = np.zeros(nv); r1[y] = 1; r1[x] = -1; r1[M + 1 + k] = -G      # t_y - t_x >= g - G(1-o)
        rows.append(r1); lo.append(g - G); hi.append(np.inf)
        r2 = np.zeros(nv); r2[x] = 1; r2[y] = -1; r2[M + 1 + k] = G       # t_x - t_y >= g - G o
        rows.append(r2); lo.append(g); hi.append(np.inf)
    for x in range(M):
        r3 = np.zeros(nv); r3[M] = 1; r3[x] = -1                          # L >= t_x
        rows.append(r3); lo.append(0); hi.append(np.inf)
    lb = np.zeros(nv); ub = np.concatenate([np.full(M + 1, np.inf), np.ones(nP)])
    lb[0] = 0; ub[0] = 0                                                  # t_0 = 0 (translation)... not WLOG for order
    ub[0] = np.inf
    res = milp(cost, constraints=[LinearConstraint(np.array(rows), np.array(lo), np.array(hi))],
               integrality=np.concatenate([np.zeros(M + 1), np.ones(nP)]), bounds=Bounds(lb, ub),
               options=dict(time_limit=time_limit, disp=False))
    if res.x is None:
        return None, None
    return res.fun, res.x[:M]
