"""§4.2 (finishing session): ILS (fls.best_ils) for W_sc on genuine D at M = 12, 20, 24, and an
actual-prime realisation with the sign matrix fixed.  Same as ../ils_run.py except:
* primes are generated on demand (smallest prime = 1 mod 4 above max(X, e^(w - 0.02)) not yet used),
  so heavy columns no longer overrun a fixed prime list (the crash of the first run);
* W_cut and its certified dual bound are read from ../perron_check.json or ../q_lean_M*.jsonl;
* every realised profile is checked exactly (verify_profile.check_pairs) and characters of support <= 4
  are evaluated in floating point (q_run_large.char_small_float).
Usage: python3 ils2_run.py M t1,t2,.. Xf kicks npool
"""
import sys, os, json, math, time
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from scipy.optimize import linprog          # noqa: E402
from fls import best_ils, sylvester_rows    # noqa: E402
from qcore import pair_index, s_to_mask     # noqa: E402
from verify_profile import check_pairs      # noqa: E402
from q_run_large import char_small_float    # noqa: E402
from gauss7 import is_prime, paley_normalised   # noqa: E402

UP = os.path.dirname(HERE)


def next_p1(n):
    """least prime p = 1 mod 4 with p >= n"""
    p = max(n, 5)
    p += (1 - p) % 4
    while not is_prime(p):
        p += 4
    return p


def lp_floor(B, S, floors):
    M, r = S.shape
    iu = pair_index(M)
    P = (S[iu[0], :] != S[iu[1], :]).astype(float) - 0.5
    res = linprog(np.ones(r), A_ub=-P, b_ub=-B[iu], bounds=[(f, None) for f in floors], method='highs')
    return (res.x, res.fun) if res.status == 0 else (None, None)


def assign_primes(w, X):
    used, out = set(), [None] * len(w)
    for j in np.argsort(w):
        p = next_p1(max(X + 1, int(math.ceil(math.exp(w[j] - 0.02)))))
        while p in used:
            p = next_p1(p + 1)
        used.add(p); out[j] = p
    return out


def realise_fixed_S(D, S, mod, X):
    """delta is added to pairs with D > 0; pairs with D = 0 get delta0 in (delta, delta/4, 0.05, 0),
    the first value for which the LP is feasible (some S admit no strictly positive margin on certain
    D = 0 pairs).  The exact integer check decides."""
    M, r = S.shape
    pos = (D > 0).astype(float)
    off = 1 - np.eye(M)
    for delta in (0.25, 0.75, 1.5, 3.0):
        floors = np.full(r, math.log(next_p1(X + 1)))
        B = None
        for d0 in (delta, delta / 4, 0.05, 0.0):
            Bt = D + delta * pos * off + d0 * (1 - pos) * off
            if lp_floor(Bt, S, floors)[0] is not None:
                B = Bt
                break
        if B is None:
            continue
        assign = None
        for _ in range(30):
            w, val = lp_floor(B, S, floors)
            if w is None:
                break
            new = assign_primes(w, X)
            if new == assign:
                break
            assign = new
            floors = np.array([math.log(p) - 0.02 for p in assign])
        if assign is None:
            continue
        cols = [(s_to_mask(S[:, j]), [assign[j]]) for j in range(r)]
        ok, worst = check_pairs(cols, mod, X)
        if ok:
            return cols, delta, worst
    return None, None, None


def stair(X, r):
    out, p = [], X + 1
    while len(out) < r:
        p = next_p1(p); out.append(p); p += 1
    return float(sum(math.log(q) for q in out))


def wcut_lookup(M, t, X):
    for rec in json.load(open(os.path.join(UP, 'perron_check.json'))):
        if rec['M'] == M and rec['t'] == t and rec['X'] == X:
            return rec['W'], rec['lb']
    for fn in ('q_lean_M%d.jsonl' % M,):
        p = os.path.join(UP, fn)
        if os.path.exists(p):
            for line in open(p):
                rec = json.loads(line)
                if rec['t'] == t and rec['X'] == X:
                    return rec['Wcut'], rec['Wcut_lb']
    return None, None


def run(args):
    M, t, Xf, kicks = args
    t0 = time.time()
    D = np.load(os.path.join(UP, 'data/D_construct_M%d_t%d_u1_X%d.npy' % (M, t, Xf)))
    mod = np.load(os.path.join(UP, 'data/mod_construct_M%d_t%d_u1_X%d.npy' % (M, t, Xf)))
    X = Xf * M
    w0 = math.log(next_p1(X + 1))
    H = sylvester_rows(M)
    P = paley_normalised(M - 1)[:, 1:]
    inits = [('syl2', np.hstack([H, H])), ('syl3', np.hstack([H, H, H])),
             ('pal2', np.hstack([P, P])), ('pal3', np.hstack([P, P, P]))]
    b = best_ils(D, w0, inits, kicks=kicks, seed=t)
    wc, wlb = wcut_lookup(M, t, X)
    lam2 = 2 * float(np.linalg.eigvalsh(D)[-1])
    rec = dict(M=M, t=t, X=X, w0=w0, lam2=lam2, Wcut=wc, Wcut_lb=wlb, Dmax=float(D.max()),
               unif=2 * (M - 1) * float(D.max()), W_ils=b['W'], ils_r=int(b['r']), ils_init=b['init'],
               ils_floor_cols=int((b['w'] < w0 + 1e-6).sum()))
    cols, delta, worst = realise_fixed_S(D, b['S'], mod, X)
    if cols is not None:
        rec['prime_ils'] = sum(math.log(p) for _, pl in cols for p in pl)
        rec['prime_ils_r'] = len(cols); rec['prime_ils_delta'] = delta; rec['prime_ils_ok'] = True
        rec['prime_ils_maxp'] = max(p for _, pl in cols for p in pl)
        rec['stair'] = stair(X, len(cols))
        cnt, cmin = char_small_float(cols, M, X)
        rec['char4_count'] = cnt; rec['char4_min'] = cmin
        rec['prime_ils_cols'] = [[int(m), [int(p) for p in pl]] for m, pl in cols]
    rec['S'] = b['S'].astype(int).tolist(); rec['w'] = [float(v) for v in b['w']]
    rec['time'] = time.time() - t0
    return rec


if __name__ == '__main__':
    M = int(sys.argv[1]); trials = [int(v) for v in sys.argv[2].split(',')]
    Xfs = [int(v) for v in sys.argv[3].split(',')]; kicks = int(sys.argv[4]); npool = int(sys.argv[5])
    jobs = [(M, t, Xf, kicks) for t in trials for Xf in Xfs]
    with Pool(min(npool, len(jobs))) as pool, open(os.path.join(HERE, 'ils2_M%d.jsonl' % M), 'a') as f:
        for rec in pool.imap_unordered(run, jobs):
            f.write(json.dumps(rec) + '\n'); f.flush()
            print('M=%d t=%d X=%d: Wcut %s W_ils %.2f (r %d, %d at floor, %s) prime %s (delta %s) %.0fs' % (
                rec['M'], rec['t'], rec['X'], rec['Wcut'], rec['W_ils'], rec['ils_r'], rec['ils_floor_cols'],
                rec['ils_init'], rec.get('prime_ils'), rec.get('prime_ils_delta'), rec['time']), flush=True)
