"""MILP with lazy constraints for the normalised problem

  Omega(M; rho0, rho1) = min sum_S om_S
     over om in ({0} cup [1, U])^{patterns}, windings kappa in Z^{M-1},
     s.t. sigma_lambda(om) >= r(lambda) = rho0 + rho1 ||lambda||_1/2 for every relevant lambda
          (lambda in L_T \\ {0}, T = supp om, lambda.(0,kappa) in Z).

Formulation (all constraints valid for every feasible point with Omega <= Ubar, so the optimum over
any finite pool of them is a lower bound for min(Omega, Ubar)):
  * integral c:   sum_S om_S (|c(S)| - 1/2) >= r(c);
  * rational class C (lambda mod Z^M, denominator D, numerator vector l = D*frac):
      l.kappa - D q_C - s_C = 0,  0 <= s_C <= D-1 integer,  z_C binary,  z_C <= s_C;
      for each stored lift lambda of C:
        sum_{S good} om_S (|lambda(S)| - 1/2) + Bm sum_{S bad} om_S + Bm z_C >= r(lambda),
      S 'bad' iff lambda(S) is not an integer (then S in T makes lambda irrelevant; om_S >= 1),
      Bm = r(lambda) + Ubar/2.  z_C = 1 is possible only if lambda.kappa is not an integer.
  * sum_S om_S <= Ubar.
Separation is exact (verify.py) on the rationalised solution.  Termination: the rationalised MILP
solution passes the exact verifier; its value is then optimal up to the MILP solver's tolerances
(HiGHS, double precision): the reported lower bound is the solver's dual bound on the final pool.

usage: python3 milp.py M rho0 rho1 Ubar [time_limit] [strict]
"""
import sys, time, math, json, itertools
from fractions import Fraction
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix, vstack as svstack
from core import all_patterns, Qmat, rows_of, lcm_list, frac_part
from verify import verify


def ternary(M, max_supp=None):
    out = []
    for c in itertools.product((-1, 0, 1), repeat=M):
        if sum(c) != 0 or not any(c):
            continue
        if next(v for v in c if v) < 0:
            continue
        if max_supp and sum(1 for v in c if v) > max_supp:
            continue
        out.append(c)
    return np.array(out, dtype=np.int64)


def multiset_perms(seq):
    items = sorted(seq)
    n = len(items)
    used = [False] * n
    cur = []
    def rec():
        if len(cur) == n:
            yield tuple(cur)
            return
        prev = None
        for i in range(n):
            if used[i] or (prev is not None and items[i] == prev):
                continue
            used[i] = True
            cur.append(items[i])
            yield from rec()
            cur.pop()
            used[i] = False
            prev = items[i]
    yield from rec()


def orbit(vec, cap=3000, seed=0):
    """Distinct images of vec under all permutations of coordinates (all of them if there are at most
    cap, else cap random ones plus vec itself)."""
    from collections import Counter
    import random
    cnt = Counter(vec)
    tot = math.factorial(len(vec))
    for v in cnt.values():
        tot //= math.factorial(v)
    if tot <= cap:
        return list(multiset_perms(vec))
    rng = random.Random(seed)
    out = {tuple(vec)}
    v = list(vec)
    for _ in range(cap):
        rng.shuffle(v)
        out.add(tuple(v))
    return list(out)


class Problem:
    def __init__(self, M, rho0, rho1, Ubar, strict=False, init_supp=6):
        self.M = M
        self.rho0 = Fraction(rho0)
        self.rho1 = Fraction(rho1)
        self.Ubar = float(Ubar)
        self.strict = strict
        self.pats = all_patterns(M)
        self.n = len(self.pats)
        self.Qall = Qmat(M, self.pats)                     # n x (M-1)
        C = ternary(M, init_supp)
        self.int_rows = []                                  # list of (coef row, rhs)
        self.cover_rows = []                                # v(c) != 0: sum_{S: c(S) != 0} om_S >= 1
        self.int_keys = set()
        for c in C:
            self.add_integral(c)
        self.classes = {}                                   # key -> dict
        self.order = []

    # r(lambda) for lambda given as full vector (length M) of Fractions or ints
    def req(self, lam_full):
        return float(self.rho0 + self.rho1 * sum(abs(Fraction(v)) for v in lam_full) / 2)

    def add_integral(self, c):
        c = np.array([int(v) for v in c], dtype=np.int64)
        if c[np.nonzero(c)[0][0]] < 0:
            c = -c
        key = tuple(c)
        if key in self.int_keys:
            return False
        self.int_keys.add(key)
        mu = c[1:]
        vals = np.abs(self.Qall @ mu).astype(float) - 0.5
        self.int_rows.append((vals, self.req(c)))
        self.cover_rows.append((np.abs(self.Qall @ mu) > 0).astype(float))
        return True

    def add_rational(self, mu_frac):
        """mu_frac: list of Fractions (lambda_1..lambda_{M-1}); adds its class and this lift."""
        f = frac_part(mu_frac)
        D = lcm_list([v.denominator for v in f])
        if D == 1:
            lam = [-sum(mu_frac)] + list(mu_frac)
            return self.add_integral([int(v) for v in lam])
        key = tuple(f)
        if key not in self.classes:
            l = [int(v * D) for v in f]
            self.classes[key] = {'D': D, 'l': l, 'lifts': [], 'liftkeys': set()}
            self.order.append(key)
        cl = self.classes[key]
        lk = tuple(mu_frac)
        if lk in cl['liftkeys']:
            return False
        cl['liftkeys'].add(lk)
        Dm = np.array([int(v * D) for v in mu_frac], dtype=np.int64)
        vals = self.Qall @ Dm                                # D lambda(S)
        good = (vals % D == 0)
        lam = [-sum(mu_frac)] + list(mu_frac)
        r = self.req(lam)
        Bm = r + self.Ubar / 2 + 1.0
        row = np.where(good, np.abs(vals) / D - 0.5, Bm)
        cl['lifts'].append((row, r, Bm))
        return True

    def solve(self, time_limit=600, verbose=False):
        n = self.n
        M = self.M
        nk = M - 1
        ncl = len(self.order)
        L = lcm_list([self.classes[k]['D'] for k in self.order]) if ncl else 1
        nv = n + nk + 3 * ncl
        cons = []
        # integral rows
        A = np.array([r for r, _ in self.int_rows])
        b = np.array([rhs for _, rhs in self.int_rows])
        if self.strict:
            b = b + 1e-6
        A1 = np.zeros((A.shape[0], nv)); A1[:, :n] = A
        cons.append(LinearConstraint(A1, lb=b, ub=np.inf))
        Cv = np.array(self.cover_rows)
        A2 = np.zeros((Cv.shape[0], nv)); A2[:, :n] = Cv
        cons.append(LinearConstraint(A2, lb=1.0, ub=np.inf))
        a = np.zeros(nv); a[:n] = 1
        cons.append(LinearConstraint(a[None, :], lb=-np.inf, ub=self.Ubar))
        lo = np.zeros(nv); hi = np.zeros(nv); integ = np.zeros(nv)
        lo[:n] = 1; hi[:n] = self.Ubar; integ[:n] = 2               # semi-continuous
        lo[n:n + nk] = 0; hi[n:n + nk] = max(L - 1, 0); integ[n:n + nk] = 1
        if ncl:
            nlift = sum(len(self.classes[k]['lifts']) for k in self.order)
            Ar = lil_matrix((nlift, nv)); br = np.zeros(nlift)
            Am = lil_matrix((2 * ncl, nv))
            i = 0
            for c, key in enumerate(self.order):
                cl = self.classes[key]
                base = n + nk + 3 * c
                for row, r, Bm in cl['lifts']:
                    nz = np.nonzero(row)[0]
                    for j in nz:
                        Ar[i, j] = row[j]
                    Ar[i, base + 2] = Bm
                    br[i] = r + (1e-6 if self.strict else 0)
                    i += 1
                D = cl['D']
                for x in range(nk):
                    if cl['l'][x]:
                        Am[2 * c, n + x] = cl['l'][x]
                Am[2 * c, base] = -D
                Am[2 * c, base + 1] = -1
                Am[2 * c + 1, base + 2] = 1
                Am[2 * c + 1, base + 1] = -1
                integ[base] = 1; lo[base] = -(sum(cl['l']) * L) // D - 2; hi[base] = (sum(cl['l']) * L) // D + 2
                integ[base + 1] = 1; lo[base + 1] = 0; hi[base + 1] = D - 1
                integ[base + 2] = 1; lo[base + 2] = 0; hi[base + 2] = 1
            cons.append(LinearConstraint(csr_matrix(Ar), lb=br, ub=np.inf))
            lbm = np.zeros(2 * ncl); ubm = np.zeros(2 * ncl); lbm[1::2] = -np.inf
            cons.append(LinearConstraint(csr_matrix(Am), lb=lbm, ub=ubm))
        cobj = np.zeros(nv); cobj[:n] = 1
        res = milp(cobj, constraints=cons, integrality=integ, bounds=Bounds(lo, hi),
                   options={'time_limit': time_limit, 'mip_rel_gap': 1e-9, 'disp': verbose,
                            'presolve': True})
        return res, L

    def extract(self, x):
        n = self.n
        om = x[:n]
        S = [i for i in range(n) if om[i] > 0.5]
        T = [self.pats[i] for i in S]
        w = [Fraction(float(om[i])).limit_denominator(240) for i in S]
        kap = [int(round(v)) for v in x[n:n + self.M - 1]]
        return T, w, kap


def run(M, rho0, rho1, Ubar, time_limit=600, strict=False, max_iter=500, verbose=True, init_supp=6,
        log=None, orbit_cap=3000, lifts_per_class=3):
    P = Problem(M, rho0, rho1, Ubar, strict, init_supp)
    t0 = time.time()
    hist = []
    for it in range(1, max_iter + 1):
        res, L = P.solve(time_limit)
        if res.x is None:
            msg = f'MILP returned no solution: status {res.status} {res.message}'
            if verbose:
                print(msg, flush=True)
            return {'status': 'infeasible_or_error', 'message': msg, 'iter': it,
                    'dual_bound': getattr(res, 'mip_dual_bound', None)}
        T, w, kap = P.extract(res.x)
        val = float(res.fun)
        dual = getattr(res, 'mip_dual_bound', None)
        r = verify(M, T, w, kap, P.rho0, P.rho1, strict=strict, max_viol=400)
        nviol = len(r.get('violators', []))
        if verbose:
            print(f'[{time.time()-t0:7.1f}s] it {it}: obj {val:.6f} dual {dual} st {res.status} '
                  f'|T|={len(T)} G={r.get("group")} L={L} int {len(P.int_rows)} cls {len(P.order)} '
                  f'viol {nviol} min_slack {r.get("min_slack")} ok {r.get("ok")}', flush=True)
        hist.append((it, val, dual, nviol))
        if r.get('ok') is None:
            return {'status': 'verifier_inconclusive', 'reason': r.get('reason'), 'value': val}
        if r['ok']:
            status = 'optimal' if res.status == 0 else 'time_limit_feasible'
            out = {'status': status, 'M': M, 'rho0': str(P.rho0), 'rho1': str(P.rho1), 'strict': strict,
                   'value_float': val, 'Omega_exact': str(sum(w)), 'dual_bound': dual, 'iter': it,
                   'support': [rows_of(s, M) for s in T], 'weights': [str(x) for x in w], 'kappa': kap,
                   'group': r['group'], 'min_slack': str(r['min_slack']), 'argmin': [str(v) for v in r['argmin']] if r['argmin'] else None,
                   'n_int': len(P.int_rows), 'n_classes': len(P.order),
                   'time': time.time() - t0}
            if log:
                json.dump(out, open(log, 'w'), indent=1)
            return out
        added = 0
        if r.get('reason') == 'affinely dependent rows':
            mu = r['violator']
            for img in orbit([-sum(mu)] + list(mu), orbit_cap):
                added += P.add_integral([int(v) for v in img])
        # group violators by class, keep the most violated lifts, add S_M-orbits
        bycls = {}
        for s, mu in r.get('violators', []):
            key = tuple(frac_part(mu))
            bycls.setdefault(key, []).append((s, mu))
        for key, lst in bycls.items():
            lst.sort(key=lambda t: t[0])
            for s, mu in lst[:lifts_per_class]:
                lam = [-sum(mu)] + list(mu)
                for img in orbit(lam, orbit_cap):
                    mu2 = list(img[1:])
                    if all(v.denominator == 1 for v in mu2):
                        added += P.add_integral([int(v) for v in img])
                    else:
                        added += P.add_rational(mu2)
        if added == 0:
            # violators already in pool: rounding problem; report
            if verbose:
                print('  no new constraints (rounding?)', [(str(s), [str(v) for v in mu]) for s, mu in r['violators'][:3]])
            return {'status': 'stalled', 'value': val, 'iter': it}
    return {'status': 'max_iter'}


if __name__ == '__main__':
    M = int(sys.argv[1]); rho0 = Fraction(sys.argv[2]); rho1 = Fraction(sys.argv[3]); Ubar = float(sys.argv[4])
    tl = float(sys.argv[5]) if len(sys.argv) > 5 else 600
    strict = (len(sys.argv) > 6 and sys.argv[6] == 'strict')
    log = f'results/milp_M{M}_r{sys.argv[2].replace("/", "_")}_{sys.argv[3].replace("/", "_")}{"_s" if strict else ""}.json'
    import os
    os.makedirs('results', exist_ok=True)
    out = run(M, rho0, rho1, Ubar, tl, strict, log=log)
    print('RESULT', json.dumps(out, default=str))
