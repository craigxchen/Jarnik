"""Solvers for the discrete problem (Q) of p3lower.md section 5 (round 7, p3data).

Variables: a nonnegative weight lam_s on every cut s (nonconstant +-1 vector mod sign).
Pair margins mu_xy = sum_s lam_s ([s_x != s_y] - 1/2), total W = sum_s lam_s (Lemma 1.4 of p3lower.md).

  W_cut(B)  = min W  s.t. mu >= B                              (real weights; LP)
  W_sc(B)   = min W  s.t. mu >= B, lam_s in {0} u [wmin, inf)  (prime weights > X: semi-continuous)
  plus optional character rows  sum_s lam_s (|c.s| - 1)/2 >= 0  for c in a finite list.

Certified lower bound for W_cut (Theorem A(1) of p3lower.md): for any nu >= 0 with
surplus(nu) = maxcut_nu - |nu|/2 > 0, W_cut(B) >= <nu, B>/surplus(nu).  surplus is computed by the
exact max-cut enumerator (libmaxcut.so), so the bound is exact up to floating-point rounding.
"""
import math, time
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse import csr_matrix, vstack, hstack, identity, csc_matrix
from qcore import best_cuts, mask_to_s, all_cut_masks, pair_coeffs, pair_index

BIG = 1e4


def cut_value(nu, s):
    sep = (s[:, None] != s[None, :])
    return float((nu * sep).sum() / 2)


def surplus(nu):
    """exact max_s cut_nu(s) - |nu|/2 (|nu| = sum over unordered pairs)"""
    _, vals = best_cuts(nu, 1)
    return vals[0] - nu.sum() / 4.0   # nu symmetric: sum over x<y = nu.sum()/2, half of it = /4


def char_rows(chars, S):
    """rows (|c.s| - 1)/2 for characters c (list of int vectors) and cut matrix S (cuts x M)"""
    C = np.array(chars, dtype=np.int64)
    return (np.abs(C @ S.T) - 1) / 2.0


def lp_full(B, masks=None, chars=None):
    M = B.shape[0]
    if masks is None:
        masks = all_cut_masks(M)
    P, S = pair_coeffs(masks, M)
    iu = pair_index(M)
    A = -P
    b = -B[iu]
    if chars:
        R = char_rows(chars, S)
        A = np.vstack([A, -R]); b = np.concatenate([b, np.zeros(len(chars))])
    res = linprog(np.ones(len(masks)), A_ub=A, b_ub=b, bounds=(0, None), method='highs')
    assert res.status == 0, res.message
    nu = np.zeros((M, M))
    nu[iu] = -res.ineqlin.marginals[:len(iu[0])]
    nu = nu + nu.T
    return res.fun, res.x, masks, nu


def lp_colgen(B, init_masks=(), chars=None, K=60, tol=1e-9, maxit=2000, verbose=False):
    """Column generation for W_cut(B) with exact pricing.  Returns (value, lam, masks, nu, lb).
    Character rows (if any) are kept in the master; pricing then includes their duals approximately,
    so with chars the returned lb is only the pair-dual bound (still valid: chars only raise W)."""
    M = B.shape[0]
    iu = pair_index(M)
    npair = len(iu[0])
    rng = np.random.default_rng(0)
    pool = set(init_masks)
    # seed: random bisections
    while len(pool) < 4 * M:
        perm = rng.permutation(M)
        s = np.ones(M, dtype=np.int64); s[perm[:M // 2]] = -1
        if s[0] < 0: s = -s
        pool.add(int(sum(1 << x for x in range(M) if s[x] < 0)))
    masks = sorted(pool)
    it = 0
    while True:
        it += 1
        P, S = pair_coeffs(masks, M)
        nv = len(masks)
        # variables: lam (nv), artificial a (npair)
        A = np.hstack([-P, -np.eye(npair)])
        b = -B[iu]
        if chars:
            R = char_rows(chars, S)
            A = np.vstack([A, np.hstack([-R, np.zeros((len(chars), npair))])])
            b = np.concatenate([b, np.zeros(len(chars))])
        c = np.concatenate([np.ones(nv), BIG * np.ones(npair)])
        res = linprog(c, A_ub=A, b_ub=b, bounds=(0, None), method='highs')
        assert res.status == 0, res.message
        y = -res.ineqlin.marginals
        nu = np.zeros((M, M)); nu[iu] = y[:npair]; nu = nu + nu.T
        newm, vals = best_cuts(nu, K)
        nusum = nu.sum() / 4.0
        add = [m for m, v in zip(newm, vals) if v - nusum > 1 + tol and m not in pool]
        if verbose and it % 10 == 0:
            print('  it %d cols %d obj %.6f best rc %.3e' % (it, nv, res.fun, 1 - (vals[0] - nusum)))
        if not add or it >= maxit:
            break
        for m in add:
            pool.add(m)
        masks = sorted(pool)
    lam = res.x[:nv]
    art = res.x[nv:].sum()
    assert art < 1e-7, 'artificials nonzero: %g' % art
    sp = surplus(nu)
    lb = float((nu[iu] * B[iu]).sum() / sp) if sp > 0 else 0.0
    return res.fun, lam, masks, nu, lb


def sc_milp(B, masks, wmin, U=None, chars=None, time_limit=600, extra_rows=None, x0=None, gap=0.0):
    """Semi-continuous MILP over the given cut list.  lam_s in {0} u [wmin, U].
    Returns dict(status, W, lb, lam, masks)."""
    M = B.shape[0]
    iu = pair_index(M)
    P, S = pair_coeffs(masks, M)
    n = len(masks)
    if U is None:
        U = float(np.max(B)) * 2 * M
    rows, lo, hi = [], [], []
    # pair rows: P lam >= B
    rows.append(np.hstack([P, np.zeros_like(P)])); lo.append(B[iu]); hi.append(np.full(len(iu[0]), np.inf))
    if chars:
        R = char_rows(chars, S)
        rows.append(np.hstack([R, np.zeros_like(R)])); lo.append(np.zeros(len(chars))); hi.append(np.full(len(chars), np.inf))
    I = np.eye(n)
    rows.append(np.hstack([I, -wmin * I])); lo.append(np.zeros(n)); hi.append(np.full(n, np.inf))
    rows.append(np.hstack([I, -U * I])); lo.append(np.full(n, -np.inf)); hi.append(np.zeros(n))
    A = csr_matrix(np.vstack(rows))
    cons = LinearConstraint(A, np.concatenate(lo), np.concatenate(hi))
    c = np.concatenate([np.ones(n), np.zeros(n)])
    integ = np.concatenate([np.zeros(n), np.ones(n)])
    bnds = Bounds(np.zeros(2 * n), np.concatenate([np.full(n, U), np.ones(n)]))
    t0 = time.time()
    res = milp(c, constraints=[cons], integrality=integ, bounds=bnds,
               options=dict(time_limit=time_limit, mip_rel_gap=gap, disp=False))
    out = dict(status=res.status, time=time.time() - t0, masks=masks)
    out['lb'] = getattr(res, 'mip_dual_bound', None)
    if res.x is not None:
        lam = res.x[:n].copy()
        lam[lam < 1e-9] = 0
        out['lam'] = lam
        out['W'] = float(lam.sum())
    else:
        out['lam'] = None; out['W'] = None
    return out


def lp_bounded(B, masks, lower, upper, chars=None):
    """LP over given cuts with per-cut bounds; returns (status, value, lam)."""
    M = B.shape[0]
    iu = pair_index(M)
    P, S = pair_coeffs(masks, M)
    A = -P; b = -B[iu]
    if chars:
        R = char_rows(chars, S)
        A = np.vstack([A, -R]); b = np.concatenate([b, np.zeros(len(chars))])
    res = linprog(np.ones(len(masks)), A_ub=A, b_ub=b, bounds=list(zip(lower, upper)), method='highs')
    if res.status != 0:
        return res.status, None, None
    return 0, res.fun, res.x


def dive(B, masks, wmin, chars=None, rule='largest', verbose=False):
    """LP diving for the semi-continuous problem: repeatedly fix (lam >= wmin) the free cut with the
    largest fractional weight below wmin; if infeasible, ban it instead.  Returns (W, lam) or (None, None)."""
    n = len(masks)
    lower = np.zeros(n); upper = np.full(n, np.inf)
    state = np.zeros(n, dtype=int)   # 0 free, 1 fixed (>= wmin), -1 banned
    while True:
        st, val, lam = lp_bounded(B, masks, lower, upper, chars)
        if st != 0:
            return None, None
        frac = [(lam[i], i) for i in range(n) if state[i] == 0 and 1e-9 < lam[i] < wmin - 1e-9]
        if not frac:
            lam = lam.copy(); lam[lam < 1e-9] = 0
            return val, lam
        frac.sort(reverse=(rule == 'largest'))
        i = frac[0][1]
        lower[i] = wmin; state[i] = 1
        st2, _, _ = lp_bounded(B, masks, lower, upper, chars)
        if st2 != 0:
            lower[i] = 0; upper[i] = 0; state[i] = -1


def polish(B, masks, lam, wmin, extra_masks=(), chars=None, time_limit=120):
    """sc MILP restricted to supp(lam) plus extra cuts."""
    sup = [masks[i] for i in range(len(masks)) if lam[i] > 0]
    cand = sorted(set(sup) | set(extra_masks))
    U = max(4 * float(np.max(lam)), 4 * wmin)
    return sc_milp(B, cand, wmin, U=U, chars=chars, time_limit=time_limit)


def lns(B, wmin, pool_masks, nu, init=None, iters=8, cand_size=220, time_limit=60, chars=None, seed=0, verbose=False):
    """Large-neighbourhood search for the semi-continuous problem.  pool_masks: cuts to draw from
    (e.g. LP support); nu: LP duals (cuts with small reduced cost are preferred)."""
    rng = np.random.default_rng(seed)
    M = B.shape[0]
    good, vals = best_cuts(nu, 1500)
    good = list(good)
    best = init
    if best is None:
        cand = sorted(set(pool_masks) | set(good[:cand_size]))
        best = sc_milp(B, cand, wmin, U=max(8 * wmin, 2 * float(B.max())), chars=chars, time_limit=time_limit)
    hist = [best['W']]
    for it in range(iters):
        sup = [m for m, l in zip(best['masks'], best['lam']) if l > 0]
        extra = list(rng.choice(len(good), size=min(cand_size, len(good)), replace=False))
        cand = sorted(set(sup) | set(pool_masks) | {good[i] for i in extra})
        U = max(4 * float(max(best['lam'])), 8 * wmin)
        r = sc_milp(B, cand, wmin, U=U, chars=chars, time_limit=time_limit)
        if r['W'] is not None and r['W'] < best['W'] - 1e-9:
            best = r
        hist.append(best['W'])
        if verbose:
            print('   lns it %d W %.3f' % (it, best['W']))
    best['hist'] = hist
    return best


def realize_primes(B, cand_masks, plist, time_limit=300):
    """Actual-prime MILP over a candidate cut set: every candidate cut may receive a set of distinct
    primes from plist (one physical column per prime; a cut with no prime is unused); every prime is used
    at most once; minimise sum log p subject to mu >= B.  Returns (cols, status, dual_bound)."""
    M = B.shape[0]
    iu = pair_index(M)
    T = list(cand_masks)
    P, S = pair_coeffs(T, M)
    nT, nP = len(T), len(plist)
    logp = np.log(np.array(plist, dtype=float))
    nv = nT * nP
    cost = np.tile(logp, nT)
    Apair = csr_matrix(np.kron(P, logp[None, :]))
    Aprime = csr_matrix(np.kron(np.ones((1, nT)), np.eye(nP)))
    A = vstack([Apair, Aprime]).tocsr()
    lo = np.concatenate([B[iu] + 1e-7, np.full(nP, -np.inf)])
    hi = np.concatenate([np.full(len(iu[0]), np.inf), np.ones(nP)])
    res = milp(cost, constraints=[LinearConstraint(A, lo, hi)],
               integrality=np.ones(nv), bounds=Bounds(np.zeros(nv), np.ones(nv)),
               options=dict(time_limit=time_limit, disp=False))
    if res.x is None:
        return None, res.status, None
    y = np.round(res.x).reshape(nT, nP)
    cols = [(T[s], [plist[k] for k in range(nP) if y[s, k] > 0.5]) for s in range(nT)]
    cols = [c for c in cols if c[1]]
    return cols, res.status, getattr(res, 'mip_dual_bound', None)


def heuristic_cuts(nu, nstart=300, rng=None, seeds=()):
    """Multi-start 1-flip local search for large cuts of nu.  Returns dict mask -> cut value (all local maxima found)."""
    if rng is None:
        rng = np.random.default_rng(0)
    M = nu.shape[0]
    found = {}
    starts = [np.where(rng.random(M) < 0.5, 1, -1) for _ in range(nstart)] + [np.where(np.asarray(sd) >= 0, 1, -1) for sd in seeds]
    for s in starts:
        s = s.astype(np.int64).copy()
        while True:
            g = s * (nu @ s)          # gain of flipping v: sum_y nu_vy (same) - (diff) = s_v (nu s)_v
            v = int(np.argmax(g))
            if g[v] <= 1e-12:
                break
            s[v] = -s[v]
        if len(set(s)) < 2:
            continue
        if s[0] < 0:
            s = -s
        m = int(sum(1 << x for x in range(M) if s[x] < 0))
        found[m] = float((nu * (s[:, None] != s[None, :])).sum() / 2)
    return found


def lp_colgen_heur(B, init_masks=(), K=80, tol=1e-9, maxit=400, nstart=300, verbose=False):
    """Column generation with heuristic pricing (no optimality certificate): returns an upper bound on W_cut."""
    M = B.shape[0]
    iu = pair_index(M)
    npair = len(iu[0])
    rng = np.random.default_rng(1)
    pool = set(init_masks)
    while len(pool) < 4 * M:
        perm = rng.permutation(M)
        s = np.ones(M, dtype=np.int64); s[perm[:M // 2]] = -1
        if s[0] < 0: s = -s
        pool.add(int(sum(1 << x for x in range(M) if s[x] < 0)))
    it = 0
    while True:
        it += 1
        masks = sorted(pool)
        P, S = pair_coeffs(masks, M)
        nv = len(masks)
        A = np.hstack([-P, -np.eye(npair)])
        c = np.concatenate([np.ones(nv), BIG * np.ones(npair)])
        res = linprog(c, A_ub=A, b_ub=-B[iu], bounds=(0, None), method='highs')
        assert res.status == 0, res.message
        y = -res.ineqlin.marginals
        nu = np.zeros((M, M)); nu[iu] = y; nu = nu + nu.T
        ev = np.linalg.eigh(-nu)[1]
        found = heuristic_cuts(nu, nstart=nstart, rng=rng, seeds=[ev[:, -1], ev[:, -2], ev[:, -3]])
        nusum = nu.sum() / 4.0
        add = sorted([(v, m) for m, v in found.items() if v - nusum > 1 + tol and m not in pool], reverse=True)[:K]
        if verbose and it % 10 == 0:
            print('  it %d cols %d obj %.4f add %d' % (it, nv, res.fun, len(add)))
        if not add or it >= maxit:
            break
        for v, m in add:
            pool.add(m)
    lam = res.x[:nv]
    assert res.x[nv:].sum() < 1e-7
    return res.fun, lam, masks, nu
