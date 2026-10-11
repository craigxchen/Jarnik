"""Fixed-support evaluation:  V(T; rho0, rho1) = min over windings of the LP
      min sum_S om_S  s.t. om_S >= 1 (S in T),  sigma_lambda(om) >= r(lambda) for relevant lambda,
with lazy constraint generation from the exact verifier, an exact feasibility check of the final
weights (upper bound) and an exact LP-dual lower bound for the final constraint pool.

Windings: it suffices to consider homomorphisms h with inclusion-minimal kernels (a smaller kernel
means fewer constraints).  For cyclic G_T an injective h exists and only integral characters remain.
"""
import math, itertools, json
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
from core import Group, Qmat, frac_part, lcm_list, rows_of, pattern_of, all_patterns
from verify import verify


def lam_rows(M, T, lam_list):
    Q = Qmat(M, T)
    rows = []
    for mu in lam_list:
        D = lcm_list([v.denominator for v in mu])
        v = np.array([int(x * D) for x in mu], dtype=np.int64)
        rows.append([Fraction(int(abs(t)), D) for t in Q @ v])
    return rows


def lp_on_pool(M, T, pool, rho0, rho1):
    """pool: list of mu (Fractions).  Returns (float value, float om, dual y) for
    min sum om s.t. om >= 1, sum_S om_S (|lam(S)| - 1/2) >= r(lam)."""
    rows = lam_rows(M, T, pool)
    A = np.array([[float(a) - 0.5 for a in r] for r in rows])
    rhs = np.array([float(rho0 + rho1 * (sum(abs(v) for v in mu) + abs(sum(mu))) / 2) for mu in pool])
    nT = len(T)
    res = linprog(np.ones(nT), A_ub=-A, b_ub=-rhs, bounds=[(1, None)] * nT, method='highs')
    if res.status != 0:
        return None, None, None, rows
    y = -res.ineqlin.marginals if hasattr(res, 'ineqlin') else None
    return res.fun, res.x, y, rows


def exact_dual_bound(T, pool, rows, rho0, rho1, y):
    """Exact lower bound from a rationalised dual y >= 0:
    om = 1 + om',  min |T| + sum om'  s.t. A om' >= b - A 1.  Weak duality: for y >= 0 with
    A^T y <= 1:  value >= |T| + (b - A 1).y.  We rescale y so that A^T y <= 1 holds exactly."""
    nT = len(T)
    yq = [max(Fraction(0), Fraction(float(v)).limit_denominator(10 ** 6)) for v in y]
    A = [[a - Fraction(1, 2) for a in r] for r in rows]
    b = [rho0 + rho1 * (sum(abs(v) for v in mu) + abs(sum(mu))) / 2 for mu in pool]
    col = [sum(yq[i] * A[i][j] for i in range(len(A))) for j in range(nT)]
    mx = max([c for c in col] + [Fraction(0)])
    if mx > 1:
        yq = [v / mx for v in yq]
    A1 = [sum(r) for r in A]
    return nT + sum(yq[i] * (b[i] - A1[i]) for i in range(len(A)))


def evaluate_support(M, T, rho0=0, rho1=0, max_iter=60, kernels=None, verbose=False, max_group=5000):
    rho0 = Fraction(rho0); rho1 = Fraction(rho1)
    G = Group(M, T)
    if not G.full_rank:
        return {'ok': False, 'reason': 'rank'}
    if G.order > max_group:
        return {'ok': None, 'reason': f'group order {G.order} too large'}
    if kernels is None:
        if G.order == 1:
            kers = [((), frozenset())]
        else:
            kers = G.minimal_kernels()
    else:
        kers = kernels
    best = None
    for eta, K in kers:
        kap = G.kappa_of(list(eta)) if G.invs else [0] * (M - 1)
        # initial pool: pairs and ternary 4-term
        pool = []
        for x in range(M):
            for y in range(x + 1, M):
                lam = [0] * M; lam[x] = 1; lam[y] = -1
                pool.append([Fraction(v) for v in lam[1:]])
        res = None
        for it in range(max_iter):
            val, om, y, rows = lp_on_pool(M, T, pool, rho0, rho1)
            if val is None:
                res = {'value': math.inf}
                break
            w = [Fraction(float(v)).limit_denominator(1000) for v in om]
            w = [max(v, Fraction(1)) for v in w]
            r = verify(M, T, w, kap, rho0, rho1, max_viol=300)
            if r['ok'] is None:
                res = {'value': None, 'reason': r['reason']}
                break
            if r['ok']:
                lb = exact_dual_bound(T, pool, rows, rho0, rho1, y)
                res = {'value': float(sum(w)), 'Omega_exact': sum(w), 'lower_exact': lb, 'weights': w,
                       'kappa': kap, 'eta': eta, 'min_slack': r['min_slack'], 'iters': it + 1}
                break
            new = 0
            keys = set(tuple(p) for p in pool)
            for s, mu in r['violators']:
                if tuple(mu) not in keys:
                    pool.append(mu); keys.add(tuple(mu)); new += 1
            if new == 0:
                # rounding: perturb weights up slightly and retry (feasible direction unknown) -> report
                res = {'value': None, 'reason': 'rounding stall', 'float_value': val}
                break
        if verbose:
            print('  kernel', eta, '->', res.get('value') if res else None, flush=True)
        if res and res.get('value') is not None and (best is None or res['value'] < best['value']):
            best = res
    out = {'ok': best is not None and best['value'] < math.inf, 'group': G.invs, 'n_kernels': len(kers)}
    if best:
        out.update(best)
    return out


# ------------------------------------------------------------------------------------------------
# Hadamard-type supports

def sylvester(t):
    M = 2 ** t
    H = [[(-1) ** bin(x & a).count('1') for a in range(M)] for x in range(M)]
    return H


def paley(q):
    """Normalised Paley-I Hadamard matrix of order q+1 (q = 3 mod 4 prime), first column all ones."""
    chi = lambda a: 0 if a % q == 0 else (1 if pow(a, (q - 1) // 2, q) == 1 else -1)
    M = q + 1
    # rows/cols indexed 0 = infinity, 1..q = F_q elements 0..q-1
    S = [[0] * M for _ in range(M)]
    for i in range(M):
        for j in range(M):
            if i == 0 and j == 0:
                S[i][j] = 0
            elif i == 0:
                S[i][j] = 1
            elif j == 0:
                S[i][j] = -1
            else:
                S[i][j] = chi((j - 1) - (i - 1))
    H = [[S[i][j] + (1 if i == j else 0) for j in range(M)] for i in range(M)]
    # normalise: make first column all ones by row sign changes
    for i in range(M):
        if H[i][0] == -1:
            H[i] = [-v for v in H[i]]
    for i in range(M):
        for j in range(M):
            assert sum(H[i][k] * H[j][k] for k in range(M)) == (M if i == j else 0)
    return H


def core_patterns(H):
    M = len(H)
    return [pattern_of(M, [x for x in range(M) if H[x][a] == -1]) for a in range(1, M)]


def flip_pattern(M, s, x):
    rows = set(rows_of(s, M))
    rows ^= {x}
    return pattern_of(M, sorted(rows))


if __name__ == '__main__':
    for name, H in [('Sylvester4', sylvester(2)), ('Sylvester8', sylvester(3)), ('Paley8', paley(7)),
                    ('Paley12', paley(11))]:
        M = len(H)
        T = core_patterns(H)
        G = Group(M, T)
        print(name, 'core group', G.invs, 'order', G.order)
