"""Benders/covering-cut MILP for the normalised problem (windings eliminated).

  Omega(M; rho0, rho1) = min sum_S om_S over om in ({0} cup [1, U])^{patterns} such that some
  homomorphism h of G_T (T = supp om) has sigma_lambda >= r(lambda) for all lambda in L_T \\ {0}
  with h(lambda) = 0, where r = rho0 + rho1 ||lambda||_1/2  (strict: '>').

Outer MILP (valid for every feasible point with Omega <= Ubar):
  (I)  integral c:  sum_S om_S (|c(S)| - 1/2) >= r(c),   sum_{S: c(S) != 0} om_S >= 1;
  (K)  covering cut for a covering family F = {C_1..C_m} of non-integral classes (inner.py Lemma):
         sum_i y_i + sum_{S notin Good(F)} om_S >= 1,       Good(F) = {S : C_i(S) in Z for all i},
         for every stored lift lambda of C_i:  row_lambda . om - Bm (1 - y_i) >= r(lambda),
       row_lambda(S) = |lambda(S)| - 1/2 on Good(C_i), = Bm off Good(C_i);  y_i binary.
       Validity: if T' is inside Good(F) then F lies in G_T', so every homomorphism of G_T' kills some
       C_i, which must then satisfy all its constraints.
  (S)  sum om <= Ubar.
Inner (exact): inner.analyse finds all bad classes of the rationalised MILP solution; if a
homomorphism avoids them all, the solution is feasible (windings from G.kappa_of) and optimal for the
problem up to HiGHS tolerances (the cut pool is a relaxation).  Otherwise the bad classes cover; a
minimal covering subfamily and its S_M-images are added as cuts.

usage: python3 milp2.py M rho0 rho1 Ubar [time_limit_per_solve] [strict] [img_cap]
"""
import sys, time, math, json, os, random
from fractions import Fraction
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix, csr_matrix
from core import all_patterns, Qmat, rows_of, lcm_list, frac_part
from inner import analyse, find_hom, minimal_cover
from milp import ternary, orbit


def perm_apply(lam, perm):
    out = [None] * len(lam)
    for x, v in enumerate(lam):
        out[perm[x]] = v
    return out


class Outer:
    def __init__(self, M, rho0, rho1, Ubar, strict=False, init_supp=6, eps_strict=Fraction(1, 10 ** 4)):
        self.M = M
        self.rho0 = Fraction(rho0); self.rho1 = Fraction(rho1)
        self.Ubar = float(Ubar)
        self.strict = strict
        self.eps = float(eps_strict) if strict else 0.0
        self.pats = all_patterns(M)
        self.n = len(self.pats)
        self.Qall = Qmat(M, self.pats)
        self.int_rows = []; self.cover_rows = []; self.int_keys = set()
        for c in ternary(M, init_supp):
            self.add_integral(c)
        self.lifts = {}        # class key -> list of (row, r) ; key = tuple frac (mu coords)
        self.liftkeys = {}
        self.good = {}         # class key -> bool mask over patterns
        self.cuts = []         # list of (tuple of class keys)
        self.cutkeys = set()

    def req(self, lam_full):
        return float(self.rho0 + self.rho1 * sum(abs(Fraction(v)) for v in lam_full) / 2)

    def add_integral(self, c):
        c = np.array([int(v) for v in c], dtype=np.int64)
        if c[np.nonzero(c)[0][0]] < 0:
            c = -c
        key = tuple(c)
        if key in self.int_keys:
            return 0
        self.int_keys.add(key)
        vals = np.abs(self.Qall @ c[1:])
        self.int_rows.append((vals.astype(float) - 0.5, self.req(c)))
        self.cover_rows.append((vals > 0).astype(float))
        return 1

    def add_lift(self, mu):
        f = tuple(frac_part(mu))
        D = lcm_list([v.denominator for v in f])
        if f not in self.good:
            Dm = np.array([int(v * D) for v in f], dtype=np.int64)
            self.good[f] = ((self.Qall @ Dm) % D == 0)
            self.lifts[f] = []
            self.liftkeys[f] = set()
        k = tuple(mu)
        if k in self.liftkeys[f]:
            return f, 0
        self.liftkeys[f].add(k)
        Dm = np.array([int(v * D) for v in mu], dtype=np.int64)
        vals = self.Qall @ Dm
        lam = [-sum(mu)] + list(mu)
        r = self.req(lam)
        Bm = r + self.Ubar / 2 + 1.0
        row = np.where(self.good[f], np.abs(vals) / D - 0.5, Bm)
        self.lifts[f].append((row, r, Bm))
        return f, 1

    def add_cut(self, keys):
        k = tuple(sorted(keys))
        if k in self.cutkeys:
            return 0
        self.cutkeys.add(k)
        self.cuts.append(k)
        return 1

    def solve(self, time_limit):
        n = self.n
        cut_y = []
        ny = 0
        for k in self.cuts:
            cut_y.append(list(range(ny, ny + len(k))))
            ny += len(k)
        nv = n + ny
        cons = []
        A = np.array([r for r, _ in self.int_rows]); b = np.array([rhs for _, rhs in self.int_rows]) + self.eps
        A1 = np.zeros((A.shape[0], nv)); A1[:, :n] = A
        cons.append(LinearConstraint(A1, lb=b, ub=np.inf))
        Cv = np.array(self.cover_rows)
        A2 = np.zeros((Cv.shape[0], nv)); A2[:, :n] = Cv
        cons.append(LinearConstraint(A2, lb=1.0, ub=np.inf))
        a = np.zeros(nv); a[:n] = 1
        cons.append(LinearConstraint(a[None, :], lb=-np.inf, ub=self.Ubar))
        if self.cuts:
            nrows = len(self.cuts) + sum(len(self.lifts[c]) for k in self.cuts for c in k)
            Ak = lil_matrix((nrows, nv)); bk = np.zeros(nrows)
            i = 0
            for k, ys in zip(self.cuts, cut_y):
                gmask = np.ones(n, dtype=bool)
                for c in k:
                    gmask &= self.good[c]
                for j in np.nonzero(~gmask)[0]:
                    Ak[i, j] = 1.0
                for y in ys:
                    Ak[i, n + y] = 1.0
                bk[i] = 1.0
                i += 1
                for c, y in zip(k, ys):
                    for row, r, Bm in self.lifts[c]:
                        for j in np.nonzero(row)[0]:
                            Ak[i, j] = row[j]
                        Ak[i, n + y] = -Bm
                        bk[i] = r - Bm + self.eps
                        i += 1
            cons.append(LinearConstraint(csr_matrix(Ak), lb=bk, ub=np.inf))
        lo = np.zeros(nv); hi = np.ones(nv); integ = np.ones(nv)
        lo[:n] = 1; hi[:n] = self.Ubar; integ[:n] = 2
        cobj = np.zeros(nv); cobj[:n] = 1
        res = milp(cobj, constraints=cons, integrality=integ, bounds=Bounds(lo, hi),
                   options={'time_limit': time_limit, 'mip_rel_gap': 1e-9, 'disp': False})
        return res

    def extract(self, x):
        om = x[:self.n]
        S = [i for i in range(self.n) if om[i] > 0.5]
        T = [self.pats[i] for i in S]
        w = [max(Fraction(1), Fraction(float(om[i])).limit_denominator(1000)) for i in S]
        return T, w


def run(M, rho0, rho1, Ubar, time_limit=3600, strict=False, max_iter=2000, img_cap=200, verbose=True,
        log=None, seed=0):
    P = Outer(M, rho0, rho1, Ubar, strict, init_supp=int(os.environ.get("INIT_SUPP", "6")))
    rng = random.Random(seed)
    t0 = time.time()
    for it in range(1, max_iter + 1):
        res = P.solve(time_limit)
        if res.x is None:
            out = {'status': 'no_solution', 'message': str(res.message), 'iter': it,
                   'dual_bound': getattr(res, 'mip_dual_bound', None)}
            if verbose:
                print(out, flush=True)
            if log:
                json.dump(out, open(log, 'w'), indent=1, default=str)
            return out
        T, w = P.extract(res.x)
        val = float(res.fun)
        dual = getattr(res, 'mip_dual_bound', None)
        A = analyse(M, T, w, P.rho0, P.rho1, strict)
        G = A['G']
        msg = ''
        added = 0
        if 'rank_violator' in A:
            mu = A['rank_violator']
            for img in orbit([-sum(mu)] + list(mu), img_cap):
                added += P.add_integral(img)
            msg = 'rank'
        elif A.get('too_large') or A.get('inconclusive'):
            out = {'status': 'inner_failed', 'reason': 'too_large' if A.get('too_large') else 'inconclusive'}
            print(out); return out
        else:
            for s, mu in A['int_viol']:
                lam = [-sum(mu)] + list(mu)
                for img in orbit(lam, img_cap):
                    added += P.add_integral([int(v) for v in img])
            bad = A['bad']
            eta = None
            if not A['int_viol']:
                eta = find_hom(G, list(bad.keys()))
            if eta is not None and not A['int_viol']:
                kap = G.kappa_of(eta) if G.invs else [0] * (M - 1)
                status = 'optimal' if res.status == 0 else 'feasible_time_limit'
                out = {'status': status, 'M': M, 'rho0': str(P.rho0), 'rho1': str(P.rho1), 'strict': strict,
                       'Omega_exact': str(sum(w)), 'value_float': val, 'dual_bound': dual, 'iter': it,
                       'support': [rows_of(s, M) for s in T], 'weights': [str(x) for x in w],
                       'kappa': kap, 'eta': eta, 'group': G.invs, 'n_bad_classes': len(bad),
                       'min_slack_all_classes': str(A['min_slack'][0]) if A['min_slack'] else None,
                       'n_int': len(P.int_rows), 'n_cuts': len(P.cuts), 'time': time.time() - t0}
                if verbose:
                    print(f'[{time.time()-t0:7.1f}s] it {it}: obj {val:.6f} dual {dual} |T|={len(T)} G={G.invs} FEASIBLE', flush=True)
                if log:
                    json.dump(out, open(log, 'w'), indent=1, default=str)
                return out
            if bad and not A['int_viol']:
                fam = minimal_cover(G, list(bad.keys()))
                # lifts for the family classes
                keys = []
                lams = []
                for nv in fam:
                    f, lifts = bad[nv]
                    keys.append(f)
                    lams.append([[-sum(mu)] + list(mu) for s, mu in lifts])
                # identity + random permutation images of the whole family
                perms = [list(range(M))]
                for _ in range(img_cap):
                    p = list(range(M)); rng.shuffle(p); perms.append(p)
                seen = set()
                for p in perms:
                    ks = []
                    for lift_list in lams:
                        kk = None
                        for lam in lift_list:
                            l2 = perm_apply(lam, p)
                            kk, a_ = P.add_lift([Fraction(v) for v in l2[1:]])
                            added += a_
                        ks.append(kk)
                    t = tuple(sorted(ks))
                    if t in seen:
                        continue
                    seen.add(t)
                    added += P.add_cut(ks)
                msg = f'cover |F|={len(fam)} of {len(bad)} bad'
        if verbose:
            print(f'[{time.time()-t0:7.1f}s] it {it}: obj {val:.6f} dual {dual} st {res.status} |T|={len(T)} '
                  f'G={G.invs if G.full_rank else None} int {len(P.int_rows)} cuts {len(P.cuts)} '
                  f'int_viol {len(A.get("int_viol", []))} {msg} added {added}', flush=True)
        if added == 0:
            out = {'status': 'stalled', 'iter': it, 'value': val}
            print(out); return out
    return {'status': 'max_iter'}


if __name__ == '__main__':
    M = int(sys.argv[1]); rho0 = Fraction(sys.argv[2]); rho1 = Fraction(sys.argv[3]); Ubar = float(sys.argv[4])
    tl = float(sys.argv[5]) if len(sys.argv) > 5 else 3600
    strict = len(sys.argv) > 6 and sys.argv[6] == 'strict'
    cap = int(sys.argv[7]) if len(sys.argv) > 7 else 200
    os.makedirs('results', exist_ok=True)
    tag = f'M{M}_r{sys.argv[2].replace("/", "_")}_{sys.argv[3].replace("/", "_")}{"_s" if strict else ""}'
    out = run(M, rho0, rho1, Ubar, tl, strict, img_cap=cap, log=f'results/b_{tag}.json')
    print('RESULT', json.dumps(out, default=str))
