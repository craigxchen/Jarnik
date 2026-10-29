"""Zonotope order zeta(A) of a slab box (own code).
A zonotope form of order t: c = prod_i (T_{v_i} - 1)^{a_i} delta_{k0}, sum a_i = t, P_c = Y^{k0} prod (Y^{v_i}-1)^{a_i};
it lies in V(A,t) iff all vertices k0 + sum_i s_i v_i (0 <= s_i <= a_i) lie in A (A = lattice points of the
polytope P = [-1,1]^d cap slabs, so vertices in P <=> support in A).
Classes: 'all' = primitive generators in {-2..2}^d (record/killer class; sum a_i|v_ix| <= 2 forces this range);
         'pair' = generators in L_pair = {e_x, e_x - e_y} (up to sign).
Upper bound: record LP (lpzeta.py): sum_i a_i |lam(v_i)| <= 2 b_lam for every constraint |lam| <= b_lam; a rational
dual y >= 0 with sum_lam y_lam |lam(v)| >= 1 for every candidate v gives zeta <= sum 2 b_lam y_lam (verified exactly).
Exact optimum: MILP (scipy/HiGHS) over integer x_v in {0,1,2} and k0 in {-1,0,1}^d; the witness is re-verified exactly.
"""
import itertools, math, time
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
import slib

def functionals(n, slabs):
    d = n - 1
    lams, bs = [], []
    for x in range(d):
        v = [0] * d; v[x] = 1; lams.append(v); bs.append(1)
    for S, L in slabs:
        lams.append(slib.fvec(n, S)); bs.append(L)
    return np.array(lams, dtype=np.int64), np.array(bs, dtype=np.int64)

def candidates(n, slabs, cls='all'):
    d = n - 1
    lams, bs = functionals(n, slabs)
    if cls == 'pair':
        V = np.array(slib.Lpair(d), dtype=np.int64)
    else:
        V = np.array(list(itertools.product(range(-2, 3), repeat=d)), dtype=np.int64)
        nz = np.any(V != 0, axis=1)
        V = V[nz]
        first = V[np.arange(V.shape[0]), np.argmax(V != 0, axis=1)]
        V = V[first > 0]
        g = np.gcd.reduce(np.abs(V), axis=1)
        V = V[g == 1]
    W = np.abs(V @ lams.T)
    ok = np.all(W <= 2 * bs[None, :], axis=1)
    return V[ok], lams, bs

def lp_bound(V, lams, bs):
    """exact rational upper bound for sum a_i via the dual LP (constraint generation over the candidate list)"""
    Wi = np.abs(V @ lams.T)
    caps = 2 * bs.astype(float)
    W = Wi.astype(np.float32)
    # start: candidates with small total width, plus a random sample
    rng = np.random.default_rng(0)
    act = set(np.argsort(Wi.sum(axis=1))[:2000].tolist()) | set(rng.choice(len(V), min(2000, len(V)), replace=False).tolist())
    while True:
        idx = np.array(sorted(act))
        res = linprog(caps, A_ub=-Wi[idx].astype(float), b_ub=-np.ones(len(idx)), bounds=[(0, None)] * len(caps),
                      method='highs')
        s_all = W @ res.x.astype(np.float32)
        viol = np.nonzero(s_all < 1 - 1e-7)[0]
        if len(viol) == 0:
            break
        add = viol[np.argsort(s_all[viol])[:5000]]
        act |= set(add.tolist())
    y = [Fraction(float(x)).limit_denominator(10 ** 4) for x in res.x]
    D = 1
    for f in y:
        D = D * f.denominator // math.gcd(D, f.denominator)
    if D < 10 ** 12:
        Y = np.array([int(f * D) for f in y], dtype=np.int64)
        s = Wi @ Y
    else:
        Y = np.array([int(f * D) for f in y], dtype=object)
        s = Wi.astype(object).dot(Y)
    worst = int(s.min())
    assert worst > 0
    scale = Fraction(D, worst) if worst < D else Fraction(1)
    # bound = sum caps * y * scale, verified: every v has sum_lam |lam(v)| Y_lam * scale / D >= 1
    yy = [Fraction(int(Yi), D) * scale for Yi in Y]
    bound = sum(Fraction(int(2 * b)) * t for b, t in zip(bs, yy))
    assert all(sum(Fraction(int(w)) * t for w, t in zip(row, yy)) >= 1 for row in Wi[np.argsort(s)[:50]])
    return bound, res.fun, yy

def milp_zeta(V, lams, bs, time_limit=1800):
    """maximise sum x_v s.t. lam(k0) + sum x_v pos_lam(v) <= b, -lam(k0) - sum x_v neg_lam(v) <= b, k0 in [-1,1]^d"""
    nv = V.shape[0]; d = lams.shape[1]
    Lv = V @ lams.T            # nv x m
    pos = np.maximum(Lv, 0).T  # m x nv
    neg = np.minimum(Lv, 0).T
    m = lams.shape[0]
    # variables: x (nv) then k0 (d)
    A1 = np.hstack([pos, lams])            # <= b
    A2 = np.hstack([-neg, -lams])          # <= b
    Acon = np.vstack([A1, A2]).astype(float)
    ub = np.concatenate([bs, bs]).astype(float)
    c = np.concatenate([-np.ones(nv), np.zeros(d)])
    xub = np.where(np.max(np.abs(V), axis=1) >= 2, 1, 2).astype(float)
    bounds = Bounds(np.concatenate([np.zeros(nv), -np.ones(d)]), np.concatenate([xub, np.ones(d)]))
    from scipy.sparse import csr_matrix
    res = milp(c, constraints=LinearConstraint(csr_matrix(Acon), -np.inf, ub), integrality=np.ones(nv + d),
               bounds=bounds, options={'time_limit': time_limit, 'disp': False, 'mip_rel_gap': 0})
    if res.x is None:
        return None, res
    x = np.round(res.x[:nv]).astype(int); k0 = np.round(res.x[nv:]).astype(int)
    gens = [(tuple(int(t) for t in V[i]), int(x[i])) for i in range(nv) if x[i] > 0]
    return (tuple(int(t) for t in k0), gens), res

def zono_points(k0, gens):
    """support points with coefficients of c = prod (T_v - 1)^a delta_k0 (exact integers)"""
    c = {tuple(k0): 1}
    for v, a in gens:
        for _ in range(a):
            nc = {}
            for k, val in c.items():
                k2 = tuple(ki + vi for ki, vi in zip(k, v))
                nc[k2] = nc.get(k2, 0) + val
                nc[k] = nc.get(k, 0) - val
            c = {k: val for k, val in nc.items() if val != 0}
    return c

def verify_exact(A, k0, gens, d):
    """exact checks: support in A; c in V(A,t); c not in V(A,t+1) (some order-t moment nonzero). Returns dict."""
    t = sum(a for v, a in gens)
    c = zono_points(k0, gens)
    Aset = set(tuple(int(v) for v in k) for k in A)
    # vertices: all 0/1.. combos
    verts_ok = all(k in Aset for k in c)
    # all lattice points of the zonotope (including interior multiplicity points) in A
    allpts = {tuple(k0)}
    for v, a in gens:
        allpts = {tuple(ki + s * vi for ki, vi in zip(k, v)) for k in allpts for s in range(a + 1)}
    all_ok = all(k in Aset for k in allpts)
    # moments of order < t: on the cube every monomial reduces to a standard one of no larger degree, so it
    # suffices to check all standard monomials of degree < t (exact int64 arithmetic; entries are 0/+-1).
    pts = np.array(list(c.keys()), dtype=np.int8)
    cv = np.array(list(c.values()), dtype=np.int64)
    low_ok = True
    if t > 0:
        Ev, bounds = slib.evalmat(pts, d, t - 1)
        mom = cv @ Ev.astype(np.int64)
        low_ok = bool(np.all(mom == 0))
    # order-t moment at alpha = product of the lex-leading variables of the generators (nonzero coefficient
    # of prod <v_i,u>^{a_i}); computed exactly over Z with the full monomial (not reduced)
    alpha = [0] * d
    for v, a in gens:
        alpha[next(i for i in range(d) if v[i] != 0)] += a
    s = 0
    for k, val in c.items():
        mv = 1
        for ki, ai in zip(k, alpha):
            if ai:
                mv *= ki ** ai
        s += val * mv
    top_nonzero = (s != 0)
    return dict(order=t, support_in_A=verts_ok, all_zonotope_points_in_A=all_ok, in_V_t=low_ok,
                not_in_V_t1=top_nonzero, sup=max(abs(v) for v in c.values()), npts=len(c))

# ------------------------------------------------------------------ exact DFS with LP-dual reduction
def dual_costs(V, lams, yy):
    """exact rational c_v = sum_lam y_lam |lam(v)| as (numerator array, common denominator)"""
    D = 1
    for t in yy:
        D = D * t.denominator // math.gcd(D, t.denominator)
    Y = np.array([int(t * D) for t in yy], dtype=object)
    W = np.abs(V @ lams.T).astype(object)
    num = W.dot(Y)
    return num, D

def dfs_feasible(n, slabs, V, lams, bs, yy, B, t, time_limit=600, cube_pts=None):
    """Is there a zonotope of order >= t in A?  Exact search over the candidates allowed by the dual reduction
    (sum_v x_v (c_v - 1) <= B - t).  Returns (witness or None, complete: bool, nodes)."""
    d = n - 1
    num, D = dual_costs(V, lams, yy)
    slackN = (B - t) * D            # Fraction; condition: sum x_v (num_v - D) <= slackN
    if slackN < 0:
        return None, True, 0
    exc = np.array([int(v) - D for v in num], dtype=object)
    keep = [i for i in range(len(V)) if exc[i] <= slackN]
    order = sorted(keep, key=lambda i: exc[i])
    Vk = V[order]; ek = [exc[i] for i in order]
    Lv = Vk @ lams.T                      # nk x m
    pos = np.maximum(Lv, 0); neg = np.minimum(Lv, 0)
    if cube_pts is None:
        cube_pts = slib.cube(d)
    F = cube_pts.astype(np.int64) @ lams.T  # all cube points: functional values
    base = np.all(np.abs(F) <= bs[None, :], axis=1)   # k0 itself in A
    idx0 = np.nonzero(base)[0]
    t0 = time.time(); nodes = [0]
    maxmult = np.where(np.max(np.abs(Vk), axis=1) >= 2, 1, 2)
    class Timeout(Exception): pass
    best = [None]
    def rec(start, tot, excess, mx, mn, idx, gens):
        nodes[0] += 1
        if nodes[0] % 2000 == 0 and time.time() - t0 > time_limit:
            raise Timeout
        if tot >= t:
            best[0] = (tuple(int(x) for x in cube_pts[idx[0]]), list(gens)); return True
        for j in range(start, len(Vk)):
            if excess + ek[j] > slackN:
                break                      # ek sorted ascending
            for a in range(int(maxmult[j]), 0, -1):
                ex2 = excess + a * ek[j]
                if ex2 > slackN:
                    continue
                mx2 = mx + a * pos[j]; mn2 = mn + a * neg[j]
                if np.any(mx2 - mn2 > 2 * bs):
                    continue
                Fi = F[idx]
                ok = np.all((Fi + mx2[None, :] <= bs[None, :]) & (Fi + mn2[None, :] >= -bs[None, :]), axis=1)
                idx2 = idx[ok]
                if len(idx2) == 0:
                    continue
                gens.append((tuple(int(x) for x in Vk[j]), a))
                if rec(j + 1, tot + a, ex2, mx2, mn2, idx2, gens):
                    return True
                gens.pop()
        return False
    m = lams.shape[0]
    try:
        found = rec(0, 0, 0, np.zeros(m, dtype=np.int64), np.zeros(m, dtype=np.int64), idx0, [])
    except Timeout:
        return None, False, nodes[0]
    return (best[0] if found else None), True, nodes[0]
