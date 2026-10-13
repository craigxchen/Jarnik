"""Part 2 (lean driver) for M = 12, 20, 24: the discrete problem (Q) on genuine demand matrices.

Per instance (M, trial t, X = Xf*M), genuine unit-1 demand D of a `construct` profile, C = sqrt 2:
  lam2     = 2 lambda_max(D)                         Perron lower bound
  Wcut, lb = LP optimum; lb = certified dual bound <nu,D>/surplus(nu) (exact pricing: lp_full / lp_colgen)
  unif     = 2(M-1) D_max                             uniform margins, real weights
  sc       = upper bound for the prime-floor problem (weights in {0} u [log p_min(X), inf)):
             reweighted-L1 support + restricted MILP (rwl1.py), then `iters` LNS rounds (0 in the runs
             reported in p3data.md: LNS did not improve the rwl1 solutions in tests)
  prime    = actual distinct primes = 1 mod 4, > X: floor raised to log p_k (k = n/2, n; n = #cuts of sc),
             demands shifted by delta = 0.3, 1, 2, rwl1 again, each cut rounded to the nearest unused prime
             (or product of two); then the nearest-prime menu MILP with the exact demands may lower the cost;
             every pair inequality verified in exact integers (verify_profile.check_pairs)
  stair    = sum of log p_i over the r smallest primes = 1 mod 4 above X (forced for r distinct primes)
  char4    = min over characters of support <= 4 (and 2e_x - e_y - e_z) of s_c, floating point
  Wcut_sep = LP with the separation demand for equal spacing in the greedy order (real weights only)
No separation is imposed on the discrete profiles here: Prop. 6.2 of p3data.md bounds its cost.
Usage: python3 q_lean.py M trials Xfs iters tlim npool
"""
import sys, json, math, time
import numpy as np
from multiprocessing import Pool
from qsolve import *
from qcore import primes_1mod4_above, pair_index
from q_run import sep_demand, round_nearest, profile_margins
from rwl1 import rwl1_sc, menu_milp
from q_run_large import greedy_order, positions_from_order, char_small_float
from verify_profile import check_pairs


def run_instance(args):
    M, t, Xf, iters, tlim = args
    t0 = time.time()
    D = np.load('data/D_construct_M%d_t%d_u1_X%d.npy' % (M, t, Xf))
    mod = np.load('data/mod_construct_M%d_t%d_u1_X%d.npy' % (M, t, Xf))
    X = Xf * M
    plist = primes_1mod4_above(X, 3000)
    wmin = math.log(plist[0])
    iu = pair_index(M)
    rec = dict(M=M, t=t, X=X, wmin=wmin, Dmax=float(D.max()), Dmean=float(D[iu].mean()))
    rec['lam2'] = 2 * float(np.linalg.eigvalsh(D)[-1])
    rec['unif'] = 2 * (M - 1) * float(D.max())
    if M <= 12:
        Wc, lam, masks, nu = lp_full(D)
        lb = float((nu[iu] * D[iu]).sum() / surplus(nu))
    else:
        Wc, lam, masks, nu, lb = lp_colgen(D)
    lp_sup = [masks[i] for i in range(len(masks)) if lam[i] > 1e-9]
    rec['Wcut'] = Wc; rec['Wcut_lb'] = lb; rec['lp_support'] = len(lp_sup)

    def solver(Bx):
        r0 = rwl1_sc(Bx, masks, wmin, nu=nu, time_limit=tlim, betas=(1.0, 2.0))
        return lns(Bx, wmin, lp_sup, nu, init=r0, iters=iters, time_limit=tlim) if iters > 0 else r0

    sc = solver(D)
    rec['sc'] = sc['W']; rec['sc_ncuts'] = int((sc['lam'] > 0).sum())
    # realisation with distinct actual primes: raise the floor to the k-th prime above X (so that the
    # cuts can take distinct primes near their weights), add a safety shift delta to the demands,
    # round each cut to the nearest unused prime (or product of two), check; then try to lower the
    # cost by the nearest-prime menu MILP with the exact demands.
    n0 = rec['sc_ncuts']
    cols, best_cost = None, None
    for k in (max(1, n0 // 2), n0):
        w2 = math.log(plist[k])
        for delta in (0.3, 1.0, 2.0):
            r2 = rwl1_sc(D + delta * (1 - np.eye(M)), masks, w2, nu=nu, time_limit=tlim, betas=(1.0,))
            if r2 is None or r2['W'] is None:
                continue
            c1, dev = round_nearest(r2['masks'], r2['lam'], plist)
            if not (((profile_margins(c1, M) - D)[iu] >= 1e-9).all()):
                continue
            c2 = menu_milp(D, r2['masks'], r2['lam'], plist, K=10, time_limit=120)
            for c in (c1, c2):
                if c is None:
                    continue
                cost = sum(math.log(p) for _, pl in c for p in pl)
                if best_cost is None or cost < best_cost:
                    cols, best_cost = c, cost
                    rec['realise_floor_k'] = k; rec['realise_delta'] = delta
            break
    if cols is not None:
        r_used = sum(len(pl) for _, pl in cols)
        rec['stair'] = sum(math.log(p) for p in plist[:r_used])
    if cols is not None:
        ok, worst = check_pairs(cols, mod, X)
        rec['prime'] = sum(math.log(p) for _, pl in cols for p in pl)
        rec['prime_r'] = sum(len(pl) for _, pl in cols)
        rec['prime_pairs_ok'] = ok
        cnt, cmin = char_small_float(cols, M, X)
        rec['char4_count'] = cnt; rec['char4_min'] = cmin
        rec['prime_cols'] = [[int(m), [int(p) for p in pl]] for m, pl in cols]
    tt = positions_from_order(greedy_order(D))
    Bs = sep_demand(D, tt)
    if M <= 12:
        Wcs, _, _, _ = lp_full(Bs)
    else:
        Wcs, _, _, _, _ = lp_colgen(Bs)
    rec['Wcut_sep_greedy'] = Wcs
    rec['time'] = time.time() - t0
    return rec


if __name__ == '__main__':
    M = int(sys.argv[1]); trials = [int(v) for v in sys.argv[2].split(',')]
    Xfs = [int(v) for v in sys.argv[3].split(',')]
    iters = int(sys.argv[4]); tlim = int(sys.argv[5]); npool = int(sys.argv[6])
    import os
    ov = 'q_lean_override_M%d.txt' % M          # optional "trials Xfs" override (used to shorten the M=20 run)
    if os.path.exists(ov):
        a, b = open(ov).read().split()
        trials = [int(v) for v in a.split(',')]; Xfs = [int(v) for v in b.split(',')]
    jobs = [(M, t, Xf, iters, tlim) for t in trials for Xf in Xfs]
    with Pool(min(npool, len(jobs))) as pool, open('q_lean_M%d.jsonl' % M, 'a') as f:
        for rec in pool.imap_unordered(run_instance, jobs):
            f.write(json.dumps(rec) + '\n'); f.flush()
            print('M=%d t=%d X=%d: lam2 %.1f Wcut %.1f (lb %.1f) sc %s prime %s r %s ok %s char4 %s | Wcut_sep %.1f | unif %.1f (%.0fs)' % (
                rec['M'], rec['t'], rec['X'], rec['lam2'], rec['Wcut'], rec['Wcut_lb'], rec.get('sc'), rec.get('prime'),
                rec.get('prime_r'), rec.get('prime_pairs_ok'), rec.get('char4_min'), rec['Wcut_sep_greedy'], rec['unif'], rec['time']))
            sys.stdout.flush()
