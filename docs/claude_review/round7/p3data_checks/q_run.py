"""Part 2 driver: the discrete problem (Q) on genuine demand matrices (round 7, p3data).

For an instance (M, trial t, X = Xf*M) with genuine unit-1 moduli m_xy (D = 2 log m, C = sqrt 2):
  lam2   = 2 lambda_max(D)                      Perron lower bound (Theorem A(2) of p3lower.md)
  Wcut   = LP optimum (real weights) + certified dual bound
  unif   = 2(M-1) D_max                          least W with all margins >= D_max (real weights)
  sc     = semi-continuous optimum/upper bound    weights in {0} u [log p_min(X), inf)
  prime  = actual distinct primes > X, exact verification of every pair inequality
  sep    = same with the delta-placement condition included (positions fixed, exact verification)
  unif_sep = least uniform margin with a feasible delta placement (placement MILP + bisection)
and character checks on the final profiles.  Writes one JSON line per instance to q_results.jsonl.
Usage: python3 q_run.py M [trials] [Xfs] [lns_iters] [time_per_milp]
"""
import sys, json, math, time
import numpy as np
from fractions import Fraction
from multiprocessing import Pool
from qsolve import *
from qcore import primes_1mod4_above, all_cut_masks, pair_index, margins
from verify_profile import check_pairs, check_placement, char_check_small, char_min_milp, placement_milp


def menu_realize(B, masks, lam, plist, extra=(), nopt=8, slack=0.25, time_limit=120, required=True, maxp=1, nsmall=12):
    """Each support cut chooses up to maxp primes among the nopt smallest primes with log p >= lam_s - slack
    plus the nsmall smallest primes overall (at least one prime if required); optional extra cuts may take
    up to maxp primes among the nsmall smallest primes; all primes distinct; mu >= B (+1e-7)."""
    M = B.shape[0]
    iu = pair_index(M)
    logs = np.log(np.array(plist, dtype=float))
    T = [(m, l) for m, l in zip(masks, lam) if l > 0]
    cuts, opts, req = [], [], []
    for m, l in T:
        idx = sorted(set([k for k in range(len(plist)) if logs[k] >= l - slack][:nopt]) | set(range(nsmall)))
        cuts.append(m); opts.append(idx); req.append(required)
    for m in extra:
        if m in cuts or nsmall == 0:
            continue
        cuts.append(m); opts.append(list(range(nsmall))); req.append(False)
    P, S = pair_coeffs(cuts, M)
    var = [(c, k) for c in range(len(cuts)) for k in opts[c]]
    nv = len(var)
    cost = np.array([logs[k] for c, k in var])
    Apair = np.zeros((len(iu[0]), nv))
    for v, (c, k) in enumerate(var):
        Apair[:, v] = P[:, c] * logs[k]
    rows = [Apair]; lo = [B[iu] + 1e-7]; hi = [np.full(len(iu[0]), np.inf)]
    Acut = np.zeros((len(cuts), nv))
    for v, (c, k) in enumerate(var):
        Acut[c, v] = 1
    rows.append(Acut); lo.append(np.array([1.0 if r else 0.0 for r in req])); hi.append(np.full(len(cuts), float(maxp)))
    used = sorted({k for c, k in var})
    Apr = np.zeros((len(used), nv))
    pos = {k: i for i, k in enumerate(used)}
    for v, (c, k) in enumerate(var):
        Apr[pos[k], v] = 1
    rows.append(Apr); lo.append(np.full(len(used), -np.inf)); hi.append(np.ones(len(used)))
    res = milp(cost, constraints=[LinearConstraint(csr_matrix(np.vstack(rows)), np.concatenate(lo), np.concatenate(hi))],
               integrality=np.ones(nv), bounds=Bounds(np.zeros(nv), np.ones(nv)),
               options=dict(time_limit=time_limit, disp=False))
    if res.x is None:
        return None
    cols = {}
    for v, (c, k) in enumerate(var):
        if res.x[v] > 0.5:
            cols.setdefault(cuts[c], []).append(plist[k])
    return list(cols.items())


def realize(B, sc, plist, lp_support, time_limit=120):
    attempts = ((dict(nopt=8, slack=0.25, required=True, maxp=1, nsmall=0), ()),
                (dict(nopt=10, slack=0.6, required=False, maxp=2, nsmall=12), ()),
                (dict(nopt=12, slack=1.0, required=False, maxp=3, nsmall=24), tuple(lp_support)))
    for kw, extra in attempts:
        cols = menu_realize(B, sc['masks'], sc['lam'], plist, extra=extra, time_limit=time_limit, **kw)
        if cols is not None:
            return cols
    return None


def realize_fb(B, sc, plist, lp_support, M, masks_all, nu, iters, tlim, exact):
    """realise; if impossible, raise the weight floor to the k-th prime (k = #cuts/2, #cuts, 2 #cuts) and re-solve."""
    cols = realize(B, sc, plist, lp_support)
    if cols is not None:
        return cols, math.log(plist[0])
    n = int((sc['lam'] > 0).sum())
    for k in (n // 2, n, 2 * n):
        w2 = math.log(plist[min(k, len(plist) - 1)])
        sc2 = solve_sc(B, M, w2, masks_all, lp_support, nu, max(2, iters // 2), tlim, exact)
        if sc2['W'] is None:
            continue
        cols = realize(B, sc2, plist, lp_support)
        if cols is not None:
            return cols, w2
    return None, None


def round_nearest(masks, lam, plist):
    """Give every support cut a set of distinct primes (one prime, or a product of two) whose log-sum is
    nearest to lam_s; cuts processed by decreasing weight.  Returns cols and total |deviation|."""
    import bisect
    logs = [math.log(p) for p in plist]
    used = set()
    cols, dev = [], 0.0
    order = sorted([(l, m) for m, l in zip(masks, lam) if l > 0], reverse=True)
    for l, m in order:
        best = None
        k = bisect.bisect_left(logs, l)
        # nearest unused prime on each side of the target (no fixed window: a window can be exhausted)
        lo_k = k - 1
        while lo_k >= 0 and lo_k in used:
            lo_k -= 1
        hi_k = k
        while hi_k < len(plist) and hi_k in used:
            hi_k += 1
        for kk in (lo_k, hi_k):
            if 0 <= kk < len(plist):
                d = abs(logs[kk] - l)
                if best is None or d < best[0]:
                    best = (d, (kk,))
        if l >= 2 * logs[0]:
            for i in range(len(plist)):
                if logs[i] > l - logs[0]:
                    break
                if i in used:
                    continue
                tgt = l - logs[i]
                j = bisect.bisect_left(logs, tgt)
                for jj in (j - 1, j, j + 1):
                    if 0 <= jj < len(plist) and jj != i and jj not in used:
                        d = abs(logs[i] + logs[jj] - l)
                        if best is None or d < best[0]:
                            best = (d, (i, jj))
        for kk in best[1]:
            used.add(kk)
        cols.append((m, [plist[kk] for kk in best[1]]))
        dev += best[0]
    return cols, dev


def realize_safety(B, M, wmin, masks_all, lp_sup, nu, iters, tlim, exact, plist, log=None):
    """sc with safety delta, nearest-prime rounding, exact-in-float feasibility check; delta raised on failure."""
    iu = pair_index(M)
    for delta in (0.5, 1.0, 2.0, 3.0, 5.0):
        Bd = B + delta * (1 - np.eye(M))
        scd = solve_sc(Bd, M, wmin, masks_all, lp_sup, nu, max(2, iters // 2), tlim, exact)
        if scd['W'] is None:
            continue
        cols, dev = round_nearest(scd['masks'], scd['lam'], plist)
        mu = profile_margins(cols, M)
        ok = bool(((mu - B)[iu] >= 1e-9).all())
        if log:
            log('safety %.1f: sc %.2f rounded dev %.3f ok %s' % (delta, scd['W'], dev, ok))
        if ok:
            return cols, delta
    return None, None


def profile_margins(cols, M):
    mu = np.zeros((M, M))
    for mask, pl in cols:
        s = mask_to_s(mask, M)
        w = sum(math.log(p) for p in pl)
        mu += w * ((s[:, None] != s[None, :]) - 0.5)
    np.fill_diagonal(mu, 0)
    return mu


def rational_positions(t, scale_num=10 ** 9):
    t = np.asarray(t, dtype=float)
    t = (t - t.min()) / max(t.max() - t.min(), 1e-300)
    return [Fraction(int(round(v * scale_num)), scale_num) for v in t]


def sep_demand(D, t):
    """B' = D + 2 log(1/|t_x - t_y|) + 1e-6 for fixed positions t in [0,1] (C = sqrt 2)."""
    M = D.shape[0]
    B = D.copy()
    for x in range(M):
        for y in range(M):
            if x != y:
                B[x, y] += 2 * math.log(1.0 / abs(float(t[x]) - float(t[y]))) + 1e-6
    return B


def unif_sep(D, tl=60):
    """least uniform margin mu with a feasible placement: gamma_xy = exp((D_xy - mu)/2) <= |t_x - t_y|."""
    lo, hi = float(D.max()), float(D.max()) + 4 * math.log(D.shape[0]) + 4
    gam = np.exp((D - lo) / 2); np.fill_diagonal(gam, 0)
    L, _ = placement_milp(gam, time_limit=tl)
    if L is not None and L <= 1:
        return lo
    for _ in range(14):
        mid = (lo + hi) / 2
        gam = np.exp((D - mid) / 2); np.fill_diagonal(gam, 0)
        L, _ = placement_milp(gam, time_limit=tl)
        if L is not None and L <= 1:
            hi = mid
        else:
            lo = mid
    return hi


def solve_sc(B, M, wmin, masks_all, lp_masks, nu, iters, tlim, exact):
    if exact:
        r = sc_milp(B, masks_all, wmin, U=6 * float(B.max()) + 4 * wmin, time_limit=tlim * 4)
        r['exact'] = (r['status'] == 0)
        return r
    r = lns(B, wmin, lp_masks, nu, iters=iters, time_limit=tlim)
    r['exact'] = False
    return r


def run_instance(args):
    M, t, Xf, iters, tlim = args
    t0 = time.time()
    def log(msg):
        with open('q_progress_M%d.txt' % M, 'a') as fp:
            fp.write('t=%d X=%d %.0fs %s\n' % (t, Xf * M, time.time() - t0, msg))
    D = np.load('data/D_construct_M%d_t%d_u1_X%d.npy' % (M, t, Xf))
    mod = np.load('data/mod_construct_M%d_t%d_u1_X%d.npy' % (M, t, Xf))
    X = Xf * M
    plist = primes_1mod4_above(X, 3000)
    wmin = math.log(plist[0])
    iu = pair_index(M)
    rec = dict(M=M, t=t, X=X, wmin=wmin, Dmax=float(D.max()), Dmean=float(D[iu].mean()))
    rec['lam2'] = 2 * float(np.linalg.eigvalsh(D)[-1])
    rec['unif'] = 2 * (M - 1) * float(D.max())
    # LP
    if M <= 12:
        Wc, x, masks, nu = lp_full(D)
        sp = surplus(nu); lb = float((nu[iu] * D[iu]).sum() / sp)
        masks_all = masks
    else:
        Wc, x, masks, nu, lb = lp_colgen(D)
        masks_all = None
    lp_sup = [masks[i] for i in range(len(masks)) if x[i] > 1e-9]
    rec['Wcut'] = Wc; rec['Wcut_lb'] = lb; rec['lp_support'] = len(lp_sup)
    log('lp done')
    exact = (M <= 8)
    # semi-continuous, no separation
    sc = solve_sc(D, M, wmin, masks_all, lp_sup, nu, iters, tlim, exact)
    rec['sc'] = sc['W']; rec['sc_lb'] = sc['lb'] if exact else Wc; rec['sc_exact'] = sc['exact']
    rec['sc_ncuts'] = int((sc['lam'] > 0).sum())
    log('sc done')
    if M <= 8:
        cols, rec['prime_floor'] = realize_fb(D, sc, plist, lp_sup, M, masks_all, nu, iters, tlim, exact)
    else:
        cols, rec['prime_safety'] = realize_safety(D, M, wmin, masks_all, lp_sup, nu, iters, tlim, exact, plist, log)
    if cols is not None:
        ok, worst = check_pairs(cols, mod, X)
        rec['prime'] = sum(math.log(p) for _, pl in cols for p in pl); rec['prime_r'] = sum(len(pl) for _, pl in cols)
        rec['prime_pairs_ok'] = ok
        rec['prime_cols'] = [[int(m), [int(p) for p in pl]] for m, pl in cols]
    log('realised: %s' % (cols is not None))
    # placement for the no-sep prime profile (to choose positions)
    if cols is not None:
        mu = profile_margins(cols, M)
        gam = np.exp((D - mu) / 2); np.fill_diagonal(gam, 0)
        L, tpos = placement_milp(gam, time_limit=120)
        rec['prime_span'] = L
    else:
        tpos = None
    # separation variant: positions = normalised placement of the no-sep solution, or equally spaced in its order
    best = None
    cands = []
    if tpos is not None:
        tq = rational_positions(tpos)
        order = np.argsort(tpos)
        te = [Fraction(0)] * M
        for i, x0 in enumerate(order):
            te[x0] = Fraction(i, M - 1)
        cands = [('placement', tq), ('equal', te)]
    else:
        cands = [('equal', [Fraction(i, M - 1) for i in range(M)])]
    for name, tt in cands:
        Bs = sep_demand(D, tt)
        if M <= 12:
            Wcs, xs, ms, nus = lp_full(Bs)
        else:
            Wcs, xs, ms, nus, _ = lp_colgen(Bs)
        if best is None or Wcs < best[1]:
            best = (name, Wcs, tt, Bs, xs, ms, nus)
    name, Wcs, tt, Bs, xs, ms, nus = best
    rec['sep_positions'] = name; rec['Wcut_sep'] = Wcs
    log('sep lp done')
    lp_sup_s = [ms[i] for i in range(len(ms)) if xs[i] > 1e-9]
    scs = solve_sc(Bs, M, wmin, ms if M <= 12 else None, lp_sup_s, nus, iters, tlim, exact)
    rec['sc_sep'] = scs['W']
    if M <= 8:
        cols_s, rec['prime_sep_floor'] = realize_fb(Bs, scs, plist, lp_sup_s, M, ms if M <= 12 else None, nus, iters, tlim, exact)
    else:
        cols_s, rec['prime_sep_safety'] = realize_safety(Bs, M, wmin, ms if M <= 12 else None, lp_sup_s, nus, iters, tlim, exact, plist, log)
    if cols_s is not None:
        ok1, _ = check_pairs(cols_s, mod, X)
        ok2, bad = check_placement(cols_s, mod, X, tt)
        rec['prime_sep'] = sum(math.log(p) for _, pl in cols_s for p in pl)
        rec['prime_sep_r'] = sum(len(pl) for _, pl in cols_s)
        rec['prime_sep_pairs_ok'] = ok1; rec['prime_sep_placement_ok'] = ok2
        okc, info = char_check_small(cols_s, M, X)
        rec['char_small_ok'] = okc
        rec['char_small_min'] = info[1] if okc else None
        cm, carg = char_min_milp(cols_s, M, bound=12, time_limit=60)
        rec['char_min_all'] = cm
        rec['prime_sep_cols'] = [[int(m), [int(p) for p in pl]] for m, pl in cols_s]
        rec['positions'] = [[f.numerator, f.denominator] for f in tt]
    SM = sum((M - k) * math.log((M - k) / k) for k in range(1, M) if M - k > k)
    rec['sep_lb'] = 4.0 / M * (float(D[iu].sum()) + 2 * SM)
    log('sep realised: %s' % (cols_s is not None))
    rec['unif_sep_mu'] = unif_sep(D)
    rec['unif_sep'] = 2 * (M - 1) * rec['unif_sep_mu']
    rec['time'] = time.time() - t0
    return rec


if __name__ == '__main__':
    M = int(sys.argv[1])
    ntr = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    Xfs = [int(v) for v in sys.argv[3].split(',')] if len(sys.argv) > 3 else [2, 4]
    iters = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    tlim = int(sys.argv[5]) if len(sys.argv) > 5 else 60
    trials = [int(v) for v in sys.argv[6].split(',')] if len(sys.argv) > 6 else list(range(ntr))
    jobs = [(M, t, Xf, iters, tlim) for t in trials for Xf in Xfs]
    with Pool(4) as pool, open('q_results_M%d.jsonl' % M, 'a') as f:
        for rec in pool.imap_unordered(run_instance, jobs):
            f.write(json.dumps(rec) + '\n'); f.flush()
            print('M=%d t=%d X=%d: lam2 %.1f Wcut %.1f sc %s prime %s | sep: Wcut %.1f sc %s prime %s | unif %.1f unif_sep %.1f  (%.0fs)' % (
                rec['M'], rec['t'], rec['X'], rec['lam2'], rec['Wcut'], rec.get('sc'), rec.get('prime'),
                rec['Wcut_sep'], rec.get('sc_sep'), rec.get('prime_sep'), rec['unif'], rec['unif_sep'], rec['time']))
            sys.stdout.flush()
