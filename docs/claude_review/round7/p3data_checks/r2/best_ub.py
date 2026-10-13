"""Best upper bounds for W_sc and W_prime on genuine D (finishing session, p3data.md §4.2).

Pipeline per instance (all outputs are upper bounds; nothing here is a lower bound):
 1. structured equal-weight SA: for tau = f w0, f in FGRID, and r = r_P(tau), r_P + 2, ..., start from
    k = floor(r/(M-1)) stacked copies of a randomly row-permuted Paley Hadamard design plus random
    balanced columns, SA (sacode.c) towards d_xy >= ceil(r/2 + D_xy/tau); the first success at each
    tau is refined downwards (r - 1, r - 2, ... while SA succeeds);
 2. LP polish of the weights (fls.lp_S, floor w0) for every code found; keep the best;
 3. ILS (fls.ils, 20 kicks) from the best S;
 4. actual primes: ils2_run.realise_fixed_S (on-demand primes, exact check_pairs), characters of
    support <= 4 in floating point.
Usage: python3 best_ub.py M t1,t2 Xf1,Xf2 npool
"""
import sys, os, json, math, time
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE)); sys.path.insert(0, HERE)
from eqw import thresholds, plotkin_r, random_balanced, sa_run, exact_ok, polish   # noqa: E402
from fls import ils, lp_S                                                       # noqa: E402
from gauss7 import paley_normalised                                             # noqa: E402
from ils2_run import realise_fixed_S, stair, next_p1, wcut_lookup               # noqa: E402
from q_run_large import char_small_float                                        # noqa: E402

UP = os.path.dirname(HERE)
FGRID = (1.0, 1.25, 1.5, 1.75, 2.0, 2.5)


def sa_try(D, tau, r, P, rng, iters, tries=2):
    M = D.shape[0]
    T = thresholds(D, tau, r)
    k = r // (M - 1)
    for j in range(tries):
        perm = rng.permutation(M)
        parts = [P[perm]] * k
        if r - k * (M - 1) > 0:
            parts.append(random_balanced(M, r - k * (M - 1), rng))
        S0 = np.hstack(parts)
        v, S = sa_run(T, S0, iters, T0=0.5, seed=int(rng.integers(1 << 30)))
        if v == 0 and exact_ok(S, T):
            return S
    return None


def run(args):
    M, t, Xf, iters = args
    t0 = time.time()
    D = np.load(os.path.join(UP, 'data/D_construct_M%d_t%d_u1_X%d.npy' % (M, t, Xf)))
    mod = np.load(os.path.join(UP, 'data/mod_construct_M%d_t%d_u1_X%d.npy' % (M, t, Xf)))
    X = M * Xf
    w0 = math.log(next_p1(X + 1))
    lam2 = 2 * float(np.linalg.eigvalsh(D)[-1])
    P = paley_normalised(M - 1)[:, 1:].astype(np.int8)
    rng = np.random.default_rng(1000 * M + 10 * t + Xf)
    best, per_tau = None, []
    for f in FGRID:
        tau = f * w0
        rP = plotkin_r(D, tau)
        S_ok, r_ok = None, None
        for r in range(rP, rP + 3 * M, 2):
            S = sa_try(D, tau, r, P, rng, iters)
            if S is not None:
                S_ok, r_ok = S, r
                break
        if S_ok is None:
            per_tau.append(dict(f=f, rP=rP, r=None))
            continue
        for r in range(r_ok - 1, rP - 1, -1):          # refine downwards
            S = sa_try(D, tau, r, P, rng, iters)
            if S is None:
                break
            S_ok, r_ok = S, r
        val, w = polish(D, S_ok, w0)
        per_tau.append(dict(f=f, rP=rP, r=r_ok, rtau=r_ok * tau, polish=val))
        if val is not None and (best is None or val < best[0]):
            best = (val, f, r_ok, S_ok)
    rec = dict(M=M, t=t, X=X, w0=w0, lam2=lam2, Dmax=float(D.max()), unif=2 * (M - 1) * float(D.max()),
               per_tau=per_tau)
    wc, wlb = wcut_lookup(M, t, X)
    rec['Wcut'], rec['Wcut_lb'] = wc, wlb
    if best is None:
        rec['time'] = time.time() - t0
        return rec
    rec['W_sa'], rec['sa_f'], rec['sa_r'] = best[0], best[1], best[2]
    # ILS refinement on B0 = D + 0.05 [D = 0]: its S admits strictly positive margins on D = 0 pairs,
    # which distinct primes require (two D = 0 pairs with complementary separating sets force
    # mu = 0 exactly, impossible for distinct primes by unique factorisation).
    off = 1 - np.eye(M)
    B0 = D + 0.05 * (D <= 0) * off
    out = ils(B0, best[3].astype(np.int64), w0, kicks=20, cand=40, rng=np.random.default_rng(t))
    S_sa = best[3].astype(np.int64)
    W_ub, S_ub = best[0], S_sa
    S_fin = out['S'] if out['slack'] < 1e-7 else S_sa
    v2 = lp_S(D, S_fin, w0)
    if v2[1] is not None and v2[2] < 1e-7 and v2[0] < W_ub - 1e-9:
        W_ub, S_ub = float(v2[0]), S_fin
    v, w, sl, _ = lp_S(D, S_ub, w0)
    rec['W_ub'] = float(v); rec['ub_r'] = int(S_ub.shape[1]); rec['ub_floor_cols'] = int((w < w0 + 1e-6).sum())
    rec['W_ils_B0'] = float(out['W']) if out['slack'] < 1e-7 else None
    cols, delta, worst = realise_fixed_S(D, S_fin, mod, X)
    if cols is None and S_fin is not S_sa:
        cols, delta, worst = realise_fixed_S(D, S_sa, mod, X)
    if cols is not None:
        rec['prime'] = sum(math.log(p) for _, pl in cols for p in pl)
        rec['prime_r'] = len(cols); rec['prime_delta'] = delta; rec['prime_ok'] = True
        rec['prime_maxp'] = max(p for _, pl in cols for p in pl)
        rec['stair'] = stair(X, len(cols))
        cnt, cmin = char_small_float(cols, M, X)
        rec['char4_count'] = cnt; rec['char4_min'] = cmin
        rec['prime_cols'] = [[int(m), [int(p) for p in pl]] for m, pl in cols]
    rec['S'] = S_ub.astype(int).tolist(); rec['w'] = [float(x) for x in w]
    rec['S_fin'] = np.asarray(S_fin).astype(int).tolist()
    rec['time'] = time.time() - t0
    return rec


if __name__ == '__main__':
    M = int(sys.argv[1]); ts = [int(v) for v in sys.argv[2].split(',')]
    Xfs = [int(v) for v in sys.argv[3].split(',')]; npool = int(sys.argv[4])
    iters = int(float(sys.argv[5])) if len(sys.argv) > 5 else 10_000_000
    jobs = [(M, t, Xf, iters) for t in ts for Xf in Xfs]
    with Pool(min(npool, len(jobs))) as pool, open(os.path.join(HERE, 'best_ub_M%d.jsonl' % M), 'a') as fh:
        for rec in pool.imap_unordered(run, jobs):
            fh.write(json.dumps(rec) + '\n'); fh.flush()
            print('M=%d t=%d X=%d: 2lam %.1f Wcut %s | SA %s (f %s, r %s) -> UB %s (r %s) | prime %s | %.0fs' % (
                rec['M'], rec['t'], rec['X'], rec['lam2'], rec['Wcut'], rec.get('W_sa'), rec.get('sa_f'),
                rec.get('sa_r'), rec.get('W_ub'), rec.get('ub_r'), rec.get('prime'), rec['time']), flush=True)
