"""Faster exhaustive exact certification of zero-margin lower bounds (same claim and certificates as
exhaust.py, integer/numpy rows, lazy amplitude-2 characters, multiprocessing).

Claim certified:  no (T, om, h) with om_S >= 1 and sum om < v is admissible
   (mode 'closed': sigma >= 0 on relevant lambda; mode 'strict': sigma > 0).
For each S_M-orbit representative T with M-1 <= |T| <= ceil(v)-1 (orderly generation, exhaust.py),
(AI) is checked exactly (core.Group); for each inclusion-minimal kernel a finite family of relevant
characters is used:  all integral ternary c; all lifts with ||lambda||_inf <= 3/2 of the kernel's
rational classes; and, only if needed, all integral c with entries in [-2, 2].  Constraints are
written as integer rows  A om >= 0  (A = 2D|lambda(S)| - D).  Certificates (exact integers):
   Farkas (closed infeasible):  y >= 0,  A^T y <= 0,  y.(A 1) < 0;
   dual bound:                  y >= 0,  A^T y <= 1  =>  min sum om >= |T| - y.(A 1)  >= v;
   Gordan (strictly infeasible): y >= 0, y != 0, A^T y <= 0.
usage: python3 exhaust2.py M v mode [workers]
"""
import sys, itertools, time, math
from fractions import Fraction
from multiprocessing import Pool
import numpy as np
from scipy.optimize import linprog
from core import all_patterns, Group, frac_part, lcm_list, rows_of
from exhaust import perm_tables, canonical_supports

G_M = None; G_C3 = None; G_C5 = None; G_V = None; G_MODE = None; G_PATS = None


def init(M, v, mode):
    global G_M, G_C3, G_C5, G_V, G_MODE, G_PATS
    G_M = M; G_V = Fraction(v); G_MODE = mode
    G_PATS = all_patterns(M)
    c3, c5 = [], []
    for c in itertools.product(range(-2, 3), repeat=M):
        if sum(c) != 0 or not any(c):
            continue
        if next(a for a in c if a) < 0:
            continue
        (c3 if max(abs(a) for a in c) == 1 else c5).append(c)
    G_C3 = np.array(c3, dtype=np.int64)
    G_C5 = np.array(c5, dtype=np.int64)


def pmat(T):
    return np.array([[(s >> x) & 1 for x in range(G_M)] for s in T], dtype=np.int64)


def rat_rows(G, ker, P):
    rows = []
    for nv in ker:
        mu = G.element(nv)
        f = frac_part(mu)
        D = lcm_list([q.denominator for q in f])
        cands = []
        for x in range(G_M - 1):
            cands.append([int((f[x] + j) * D) for j in range(-2, 2) if abs(f[x] + j) <= Fraction(3, 2)])
        X = np.array(list(itertools.product(*cands)), dtype=np.int64).reshape(-1, G_M - 1)
        lam0 = -X.sum(axis=1)
        X = X[np.abs(lam0) * 2 <= 3 * D]
        lam0 = -X.sum(axis=1)
        full = np.hstack([lam0[:, None], X])
        rows.append(2 * np.abs(full @ P.T) - D)
    return np.vstack(rows) if rows else np.zeros((0, P.shape[0]), dtype=np.int64)


def ratvec(y, den=10 ** 6):
    q = [max(Fraction(0), Fraction(float(t)).limit_denominator(den)) for t in y]
    L = lcm_list([t.denominator for t in q]) if q else 1
    return np.array([int(t * L) for t in q], dtype=object), L


def certify(A, nT):
    """A: integer rows (constraints A om >= 0, om >= 1).  Returns ('infeasible'|'bound'|'gordan'|None, info)."""
    Af = A.astype(float)
    A1 = A.sum(axis=1)
    res = linprog(np.ones(nT), A_ub=-Af, b_ub=A1.astype(float), bounds=[(0, None)] * nT, method='highs')
    if res.status == 2:
        m = A.shape[0]
        r2 = linprog(np.zeros(m), A_ub=np.vstack([Af.T, A1[None, :].astype(float)]),
                     b_ub=np.concatenate([np.zeros(nT), [-1.0]]), bounds=[(0, None)] * m, method='highs')
        if r2.status == 0:
            yi, L = ratvec(r2.x)
            AT = A.T.astype(object) @ yi
            if all(t <= 0 for t in AT) and (A1.astype(object) @ yi) < 0:
                return 'infeasible', None
        return None, 'farkas failed'
    if res.status != 0:
        return None, 'lp status %d' % res.status
    y = -res.ineqlin.marginals
    yi, L = ratvec(y)
    AT = A.T.astype(object) @ yi                     # = L * (A^T y)
    mx = max([Fraction(int(t), L) for t in AT] + [Fraction(0)])
    scale = Fraction(1) if mx <= 1 else 1 / mx
    lb = nT - scale * Fraction(int(A1.astype(object) @ yi), L)
    return 'bound', lb


def gordan(A, nT):
    m = A.shape[0]
    res = linprog(np.zeros(m), A_ub=A.T.astype(float), b_ub=np.zeros(nT), A_eq=np.ones((1, m)), b_eq=[1.0],
                  bounds=[(0, None)] * m, method='highs')
    if res.status != 0:
        return False
    for den in (10 ** 4, 10 ** 6, 10 ** 9):
        yi, L = ratvec(res.x, den)
        if not any(yi):
            continue
        AT = A.T.astype(object) @ yi
        if all(t <= 0 for t in AT):
            return True
    return False


def work(sup):
    T = [G_PATS[i] for i in sup]
    G = Group(G_M, T)
    if not G.full_rank:
        return ('rank', None)
    P = pmat(T)
    A3 = 2 * np.abs(G_C3 @ P.T) - 1
    kers = [((), [])] if not G.invs else [(eta, sorted(K)) for eta, K in G.minimal_kernels()]
    out = []
    for eta, ker in kers:
        Ar = rat_rows(G, ker, P)
        A = np.vstack([A3, Ar]) if Ar.shape[0] else A3
        kind, info = certify(A, len(T))
        if kind == 'bound' and info < G_V:
            A5 = 2 * np.abs(G_C5 @ P.T) - 1
            A = np.vstack([A, A5])
            kind, info = certify(A, len(T))
        if kind == 'infeasible' or (kind == 'bound' and info >= G_V):
            out.append(kind)
            continue
        if G_MODE == 'strict' and gordan(A, len(T)):
            out.append('gordan')
            continue
        return ('FAIL', ([rows_of(s, G_M) for s in T], eta, kind, str(info)))
    return ('ok', out)


def main(M, v, mode, workers):
    pats, tabs = perm_tables(M)
    tmax = math.ceil(Fraction(v)) - 1
    t0 = time.time()
    stats = {}
    fails = []
    with Pool(workers, initializer=init, initargs=(M, v, mode)) as pool:
        for t in range(M - 1, tmax + 1):
            sups = canonical_supports(M, t, tabs)
            print(f'[{time.time()-t0:.0f}s] |T|={t}: {len(sups)} orbit representatives', flush=True)
            for i, r in enumerate(pool.imap_unordered(work, sups, chunksize=64)):
                k = r[0]
                stats[k] = stats.get(k, 0) + 1
                if k == 'ok':
                    for kk in r[1]:
                        stats['case_' + kk] = stats.get('case_' + kk, 0) + 1
                if k == 'FAIL':
                    fails.append(r[1])
                if (i + 1) % 20000 == 0:
                    print(f'   ... {i+1} done, stats {stats}', flush=True)
    print('stats', stats, f'time {time.time()-t0:.0f}s')
    print('uncertified:', len(fails))
    for f in fails[:20]:
        print('  ', f)
    print('CERTIFIED' if not fails else 'NOT CERTIFIED', f'M={M} mode={mode}: no admissible point with Omega <', v)


if __name__ == '__main__':
    main(int(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 3)
