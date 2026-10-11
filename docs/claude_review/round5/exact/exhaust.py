"""Exhaustive exact certification of lower bounds for the zero-margin problems (small M).

Claim certified:  every (T, om, h) with om_S >= 1 and sum om < v is infeasible
   (mode 'closed': sigma_lambda >= 0 for relevant lambda;  mode 'strict': sigma_lambda > 0).
Since om_S >= 1, such a support has M-1 <= |T| <= ceil(v) - 1.  Supports are enumerated up to the
action of S_M on rows (canonical form = lexicographically least sorted image).  For each support T
and each homomorphism class (inclusion-minimal kernels suffice: a larger kernel only adds
constraints), a finite family F of relevant characters is used:
   all integral ternary c, all integral c with entries in {-2..2}, and all lifts with ||lambda||_inf<=3/2
   of the kernel's rational classes.
Certificates (exact Fractions, from rationalised HiGHS duals):
   closed: LP  min |T| + sum om'  s.t.  A om' >= -A 1, om' >= 0  (om = 1 + om').
           Either a Farkas vector y >= 0, A^T y <= 0, y.(-A 1) > 0  (infeasible), or a dual vector
           y >= 0, A^T y <= 1 with |T| + y.(-A 1) >= v.
   strict: strict feasibility on F <=> exists om >= 0 with A om > 0 (scale invariance).  Either a
           Gordan vector y >= 0, y != 0, A^T y <= 0 (strictly infeasible), or the closed certificate.
If some (T, h) gets neither, it is reported (it may be a genuine solution of value < v).

usage: python3 exhaust.py M v mode
"""
import sys, itertools, time, math
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
from core import all_patterns, Group, Qmat, frac_part, lcm_list, rows_of
from milp import ternary


def perm_tables(M):
    pats = all_patterns(M)
    idx = {p: i for i, p in enumerate(pats)}
    full = (1 << M) - 1
    tabs = []
    for perm in itertools.permutations(range(M)):
        t = []
        for p in pats:
            q = 0
            for x in range(M):
                if (p >> x) & 1:
                    q |= 1 << perm[x]
            if q & 1:
                q ^= full
            t.append(idx[q])
        tabs.append(t)
    return pats, np.array(tabs, dtype=np.int32)


def _canon(tabs, sub, base):
    imgs = np.sort(tabs[:, sub], axis=1)
    codes = np.zeros(imgs.shape[0], dtype=np.int64)
    for j in range(imgs.shape[1]):
        codes = codes * base + imgs[:, j]
    i = int(np.argmin(codes))
    return tuple(int(v) for v in imgs[i])


def canonical_supports(M, t, tabs, _cache={}):
    """All S_M-orbit representatives of t-subsets of patterns, by augmentation of the (t-1)-level
    representatives (every t-subset contains a (t-1)-subset, whose orbit has a representative R;
    the image of the t-subset under the map sending that subset to R contains R)."""
    n = tabs.shape[1]
    base = n + 1
    key = (M, t)
    if key in _cache:
        return _cache[key]
    if t == 0:
        res = [()]
    else:
        prev = canonical_supports(M, t - 1, tabs)
        seen = set()
        for P in prev:
            Ps = set(P)
            for s in range(n):
                if s in Ps:
                    continue
                sub = np.array(sorted(Ps | {s}), dtype=np.int32)
                seen.add(_canon(tabs, sub, base))
        res = sorted(seen)
    _cache[key] = res
    return res


def rows_matrix(M, T, lams):
    """rows |lambda(S)| - 1/2 (Fractions) for lambda given as mu-vectors."""
    Q = Qmat(M, T)
    R = []
    for mu in lams:
        D = lcm_list([v.denominator for v in mu])
        v = np.array([int(x * D) for x in mu], dtype=np.int64)
        R.append([Fraction(int(abs(s)), D) - Fraction(1, 2) for s in Q @ v])
    return R


def rat(y, den=10 ** 6):
    return [max(Fraction(0), Fraction(float(v)).limit_denominator(den)) for v in y]


def certify_closed(R, nT, v):
    """Returns ('infeasible', y) or ('bound', value) with exact verification, or (None, None)."""
    A = np.array([[float(a) for a in r] for r in R])
    A1 = A.sum(axis=1)
    # LP: min sum om'  s.t. A om' >= -A1, om' >= 0
    res = linprog(np.ones(nT), A_ub=-A, b_ub=A1, bounds=[(0, None)] * nT, method='highs')
    if res.status == 2:
        # Farkas: y >= 0, A^T y <= 0, y.(-A1) > 0
        m = len(R)
        r2 = linprog(np.zeros(m), A_ub=np.vstack([A.T, A1[None, :]]), b_ub=np.concatenate([np.zeros(nT), [-1.0]]),
                     bounds=[(0, None)] * m, method='highs')
        if r2.status != 0:
            return None, None
        y = rat(r2.x)
        AT = [sum(y[i] * R[i][j] for i in range(m)) for j in range(nT)]
        A1e = [sum(r) for r in R]
        if all(a <= 0 for a in AT) and sum(y[i] * (-A1e[i]) for i in range(m)) > 0:
            return 'infeasible', y
        return None, None
    if res.status != 0:
        return None, None
    y = -res.ineqlin.marginals
    yq = rat(y)
    m = len(R)
    col = [sum(yq[i] * R[i][j] for i in range(m)) for j in range(nT)]
    mx = max(col + [Fraction(0)])
    if mx > 1:
        yq = [q / mx for q in yq]
    A1e = [sum(r) for r in R]
    lb = nT + sum(yq[i] * (-A1e[i]) for i in range(m))
    if lb >= v:
        return 'bound', lb
    return 'low', (lb, res.fun + nT)


def certify_strict_infeasible(R, nT):
    """Gordan: y >= 0, sum y = 1, A^T y <= 0 (exact)."""
    A = np.array([[float(a) for a in r] for r in R])
    m = len(R)
    res = linprog(np.zeros(m), A_ub=A.T, b_ub=np.zeros(nT), A_eq=np.ones((1, m)), b_eq=[1.0],
                  bounds=[(0, None)] * m, method='highs')
    if res.status != 0:
        return None
    for den in (10 ** 4, 10 ** 6, 10 ** 9):
        y = rat(res.x, den)
        if sum(y) == 0:
            continue
        AT = [sum(y[i] * R[i][j] for i in range(m)) for j in range(nT)]
        if all(a <= 0 for a in AT):
            return y
    # try support-restricted exact re-solve: keep rows with y > 0 and solve exactly by rational LP is
    # not implemented; report failure
    return None


def family(M, T, G, eta, int_chars):
    lams = [list(map(Fraction, c[1:])) for c in int_chars]
    if G.invs:
        for nv in G.elements():
            if not any(nv):
                continue
            if G.hom_value(eta, nv) != 0:
                continue
            mu = G.element(nv)
            f = frac_part(mu)
            # lifts with ||lambda||_inf <= 3/2
            cands = []
            for x in range(M - 1):
                cands.append([f[x] + j for j in range(-2, 2) if abs(f[x] + j) <= Fraction(3, 2)])
            for mu2 in itertools.product(*cands):
                if abs(sum(mu2)) <= Fraction(3, 2):
                    lams.append(list(mu2))
    return lams


def main(M, v, mode):
    v = Fraction(v)
    pats, tabs = perm_tables(M)
    tern = [tuple(int(a) for a in c) for c in ternary(M)]
    amp2 = []
    for c in itertools.product(range(-2, 3), repeat=M):
        if sum(c) == 0 and any(c) and max(abs(a) for a in c) == 2:
            if next(a for a in c if a) > 0:
                amp2.append(c)
    int_chars = tern + amp2
    tmax = math.ceil(v) - 1
    stats = {'supports': 0, 'cases': 0, 'rank': 0, 'infeasible': 0, 'bound': 0, 'strict_infeasible': 0}
    problems = []
    t0 = time.time()
    for t in range(M - 1, tmax + 1):
        sups = canonical_supports(M, t, tabs)
        print(f'|T|={t}: {len(sups)} supports up to S_{M}', flush=True)
        for sup in sups:
            T = [pats[i] for i in sup]
            stats['supports'] += 1
            G = Group(M, T)
            if not G.full_rank:
                stats['rank'] += 1
                continue
            kers = [((), None)] if not G.invs else G.minimal_kernels()
            for eta, K in kers:
                stats['cases'] += 1
                lams = family(M, T, G, list(eta), int_chars)
                R = rows_matrix(M, T, lams)
                kind, data = certify_closed(R, len(T), v)
                if kind == 'infeasible':
                    stats['infeasible'] += 1
                    continue
                if kind == 'bound':
                    stats['bound'] += 1
                    continue
                if mode == 'strict':
                    y = certify_strict_infeasible(R, len(T))
                    if y is not None:
                        stats['strict_infeasible'] += 1
                        continue
                problems.append((T, eta, kind, data))
    print('stats', stats, f'time {time.time()-t0:.1f}s', flush=True)
    print('uncertified cases:', len(problems))
    for T, eta, kind, data in problems[:20]:
        print('  ', [rows_of(s, M) for s in T], eta, kind, data)
    print('CERTIFIED' if not problems else 'NOT CERTIFIED', f'M={M} mode={mode}: no feasible point with Omega <', v)


if __name__ == '__main__':
    main(int(sys.argv[1]), sys.argv[2], sys.argv[3])
