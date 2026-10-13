"""§4.2: iterated flip local search (fls.py) for the prime-floor problem W_sc on genuine D, plus an
actual-prime realisation with the sign matrix fixed, verified in exact integers.

Per instance: W_cut (from q_lean_M*.jsonl, certified), W_ils = best prime-floor value found by ILS from
Sylvester-type starts (an upper bound for W_sc), its column count, and W_prime_ils: the columns of
the ILS solution get distinct actual primes = 1 mod 4 above X as follows.  Repeat: LP over the weights
(S fixed) with demands D + delta and per-column floors f_j; give the columns, in increasing order of
weight, the smallest unused prime with log p >= w_j - 0.02; set f_j = log p_j - 0.02; until the
assignment is stable.  Then check every pair inequality exactly (check_pairs); if it fails, raise
delta (0.25, 0.75, 1.5, 3).
Usage: python3 ils_run.py M trials Xf kicks npool
"""
import sys, json, math, time
import numpy as np
from multiprocessing import Pool
from scipy.optimize import linprog
from fls import best_ils, sylvester_rows, lp_S
from qcore import primes_1mod4_above, pair_index, s_to_mask
from verify_profile import check_pairs
from q_run_large import char_small_float


def lp_floor(B, S, floors):
    M, r = S.shape
    iu = pair_index(M)
    P = (S[iu[0], :] != S[iu[1], :]).astype(float) - 0.5
    res = linprog(np.ones(r), A_ub=-P, b_ub=-B[iu], bounds=[(f, None) for f in floors], method='highs')
    return (res.x, res.fun) if res.status == 0 else (None, None)


def realise_fixed_S(D, S, plist, mod, X):
    M, r = S.shape
    logs = np.log(np.array(plist, dtype=float))
    for delta in (0.25, 0.75, 1.5, 3.0):
        B = D + delta * (1 - np.eye(M))
        floors = np.full(r, logs[0])
        assign = None
        for _ in range(30):
            w, val = lp_floor(B, S, floors)
            if w is None:
                break
            used, new = set(), [None] * r
            for j in np.argsort(w):
                k = int(np.searchsorted(logs, w[j] - 0.02))
                while k in used:
                    k += 1
                used.add(k); new[j] = k
            if max(new) >= len(plist):          # weights beyond the prime list: give up this delta
                assign = None
                break
            if new == assign:
                break
            assign = new
            floors = np.array([logs[k] - 0.02 for k in assign])
        if assign is None:
            continue
        cols = [(s_to_mask(S[:, j]), [plist[assign[j]]]) for j in range(r)]
        ok, worst = check_pairs(cols, mod, X)
        if ok:
            return cols, delta, worst
    return None, None, None


def run(args):
    M, t, Xf, kicks = args
    t0 = time.time()
    D = np.load('data/D_construct_M%d_t%d_u1_X%d.npy' % (M, t, Xf))
    mod = np.load('data/mod_construct_M%d_t%d_u1_X%d.npy' % (M, t, Xf))
    X = Xf * M
    plist = primes_1mod4_above(X, 3000)
    w0 = math.log(plist[0])
    H = sylvester_rows(M)
    inits = [('syl2', np.hstack([H, H])), ('syl3', np.hstack([H, H, H]))]
    b = best_ils(D, w0, inits, kicks=kicks, seed=t)
    rec = dict(M=M, t=t, X=X, w0=w0, W_ils=b['W'], ils_r=int(b['r']), ils_init=b['init'],
               ils_floor_cols=int((b['w'] < w0 + 1e-6).sum()))
    cols, delta, worst = realise_fixed_S(D, b['S'], plist, mod, X)
    if cols is not None:
        rec['prime_ils'] = sum(math.log(p) for _, pl in cols for p in pl)
        rec['prime_ils_r'] = len(cols); rec['prime_ils_delta'] = delta; rec['prime_ils_ok'] = True
        rec['stair'] = float(sum(math.log(p) for p in plist[:len(cols)]))
        cnt, cmin = char_small_float(cols, M, X)
        rec['char4_count'] = cnt; rec['char4_min'] = cmin
        rec['prime_ils_cols'] = [[int(m), [int(p) for p in pl]] for m, pl in cols]
    rec['time'] = time.time() - t0
    return rec


if __name__ == '__main__':
    M = int(sys.argv[1]); trials = [int(v) for v in sys.argv[2].split(',')]
    Xfs = [int(v) for v in sys.argv[3].split(',')]; kicks = int(sys.argv[4]); npool = int(sys.argv[5])
    jobs = [(M, t, Xf, kicks) for t in trials for Xf in Xfs]
    with Pool(min(npool, len(jobs))) as pool, open('ils_M%d.jsonl' % M, 'a') as f:
        for rec in pool.imap_unordered(run, jobs):
            f.write(json.dumps(rec) + '\n'); f.flush()
            print('M=%d t=%d X=%d: W_ils %.2f (r %d, %d at floor) prime %s (delta %s) %.0fs' % (
                rec['M'], rec['t'], rec['X'], rec['W_ils'], rec['ils_r'], rec['ils_floor_cols'],
                rec.get('prime_ils'), rec.get('prime_ils_delta'), rec['time']), flush=True)
