"""Part 2 at larger M (32, 44): heuristic upper bounds for (Q); certified lower bound 2 lambda_max(D).

Same quantities as q_run.py, with these changes (all upper bounds are still verified exactly):
  * LP: column generation with heuristic pricing (lp_colgen_heur) -> an upper bound on W_cut;
    the certified lower bound is the Perron bound 2 lambda_max(D) (and Prop. 3.2 for +sep);
  * sc: large-neighbourhood search over restricted MILPs;
  * +sep: fixed positions t_i = i/(M-1) along an ordering (identity, or the order that puts the
    heaviest pairs far apart, greedy), demand D + 2 log(1/|dt|) + 1e-6, exact placement check;
  * uniform + sep: margins mu = D_max + s, smallest s on a grid of 0.25 for which equal spacing in the
    greedy order passes the exact check (an upper bound for W_unif+sep; the lower bound is W_unif);
  * characters: support <= 4 checked in floating point with exact re-check of near-zero cases; MILP min.
Usage: python3 q_run_large.py M trials Xfs iters tlim
"""
import sys, json, math, time
import numpy as np
from fractions import Fraction
from multiprocessing import Pool
from qsolve import *
from qcore import primes_1mod4_above, pair_index, mask_to_s
from q_run import realize, profile_margins, sep_demand
from verify_profile import check_pairs, check_placement, char_min_milp, basic_checks


def greedy_order(D):
    """order points so that heavy pairs are far apart: start at the endpoint of the heaviest pair, then
    repeatedly append the point with the smallest demand to the last two placed points."""
    M = D.shape[0]
    x, y = np.unravel_index(np.argmax(D), D.shape)
    order = [int(x)]
    rest = set(range(M)) - {int(x)}
    while rest:
        last = order[-2:]
        nxt = min(rest, key=lambda z: (sum(D[z, l] for l in last), z))
        order.append(nxt); rest.remove(nxt)
    return order


def positions_from_order(order):
    M = len(order)
    t = [Fraction(0)] * M
    for i, x in enumerate(order):
        t[x] = Fraction(i, M - 1)
    return t


def char_small_float(cols, M, X):
    S, ps = basic_checks(cols, M, X)
    A = (1 - S) // 2
    w = np.log(np.array(ps, dtype=float))
    W = w.sum()
    import itertools
    worst = None
    cnt = 0
    for sup in itertools.combinations(range(M), 4):
        a, b, c, d = sup
        for pos, neg in (((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))):
            v = A[pos[0]] + A[pos[1]] - A[neg[0]] - A[neg[1]]
            sc = float(np.abs(v) @ w) - W / 2
            cnt += 1
            worst = sc if worst is None else min(worst, sc)
    for x in range(M):
        for y in range(M):
            for z in range(y + 1, M):
                if x in (y, z):
                    continue
                v = 2 * A[x] - A[y] - A[z]
                sc = float(np.abs(v) @ w) - W / 2
                cnt += 1
                worst = min(worst, sc)
    return cnt, worst


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
    Wc, lam, masks, nu = lp_colgen_heur(D, maxit=60, nstart=150)
    lp_sup = [masks[i] for i in range(len(masks)) if lam[i] > 1e-9]
    rec['Wcut'] = Wc; rec['Wcut_lb'] = rec['lam2']; rec['lp_support'] = len(lp_sup)
    sc = lns(D, wmin, lp_sup, nu, iters=iters, time_limit=tlim)
    rec['sc'] = sc['W']; rec['sc_lb'] = rec['lam2']; rec['sc_exact'] = False; rec['sc_ncuts'] = int((sc['lam'] > 0).sum())
    cols = realize(D, sc, plist, lp_sup)
    if cols is not None:
        ok, _ = check_pairs(cols, mod, X)
        rec['prime'] = sum(math.log(p) for _, pl in cols for p in pl); rec['prime_r'] = sum(len(pl) for _, pl in cols)
        rec['prime_pairs_ok'] = ok
        rec['prime_cols'] = [[int(m), [int(p) for p in pl]] for m, pl in cols]
    # separation
    best = None
    for name, order in (('identity', list(range(M))), ('greedy', greedy_order(D))):
        tt = positions_from_order(order)
        Bs = sep_demand(D, tt)
        Wcs, lams, ms, nus = lp_colgen_heur(Bs, maxit=60, nstart=150)
        if best is None or Wcs < best[1]:
            best = (name, Wcs, tt, Bs, lams, ms, nus)
    name, Wcs, tt, Bs, lams, ms, nus = best
    rec['sep_positions'] = name; rec['Wcut_sep'] = Wcs
    lp_sup_s = [ms[i] for i in range(len(ms)) if lams[i] > 1e-9]
    scs = lns(Bs, wmin, lp_sup_s, nus, iters=iters, time_limit=tlim)
    rec['sc_sep'] = scs['W']
    cols_s = realize(Bs, scs, plist, lp_sup_s)
    if cols_s is not None:
        ok1, _ = check_pairs(cols_s, mod, X)
        ok2, _ = check_placement(cols_s, mod, X, tt)
        rec['prime_sep'] = sum(math.log(p) for _, pl in cols_s for p in pl); rec['prime_sep_r'] = sum(len(pl) for _, pl in cols_s)
        rec['prime_sep_pairs_ok'] = ok1; rec['prime_sep_placement_ok'] = ok2
        cnt, worst = char_small_float(cols_s, M, X)
        rec['char_small_ok'] = worst >= 1e-6; rec['char_small_min'] = worst
        cm, _ = char_min_milp(cols_s, M, bound=12, time_limit=60)
        rec['char_min_all'] = cm
        rec['prime_sep_cols'] = [[int(m), [int(p) for p in pl]] for m, pl in cols_s]
    # uniform + sep upper bound via equal spacing in greedy order
    tt = positions_from_order(greedy_order(D))
    s = 0.0
    while True:
        mu = float(D.max()) + s
        okp = True
        for x in range(M):
            for y in range(x + 1, M):
                if math.exp((D[x, y] - mu) / 2) * (1 + 1e-9) > abs(float(tt[x] - tt[y])):
                    okp = False; break
            if not okp:
                break
        if okp:
            break
        s += 0.25
    rec['unif_sep_mu'] = float(D.max()) + s
    rec['unif_sep'] = 2 * (M - 1) * rec['unif_sep_mu']
    rec['unif_sep_is_ub'] = True
    rec['time'] = time.time() - t0
    return rec


if __name__ == '__main__':
    M = int(sys.argv[1]); ntr = int(sys.argv[2])
    Xfs = [int(v) for v in sys.argv[3].split(',')]
    iters = int(sys.argv[4]); tlim = int(sys.argv[5])
    jobs = [(M, t, Xf, iters, tlim) for t in range(ntr) for Xf in Xfs]
    with Pool(min(4, len(jobs))) as pool, open('q_results_M%d.jsonl' % M, 'a') as f:
        for rec in pool.imap_unordered(run_instance, jobs):
            f.write(json.dumps(rec) + '\n'); f.flush()
            print('M=%d t=%d X=%d: lam2 %.1f Wcut<=%.1f sc %s prime %s | sep: Wcut<=%.1f prime %s | unif %.1f unif_sep<=%.1f (%.0fs)' % (
                rec['M'], rec['t'], rec['X'], rec['lam2'], rec['Wcut'], rec.get('sc'), rec.get('prime'), rec['Wcut_sep'],
                rec.get('prime_sep'), rec['unif'], rec['unif_sep'], rec['time']))
            sys.stdout.flush()
