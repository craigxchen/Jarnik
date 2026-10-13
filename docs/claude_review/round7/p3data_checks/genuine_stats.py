"""Part 1: genuine pair demands for the fully flipped Paley profiles with actual Gaussian primes.

For M in {12,20,24,32,44,60,72,84} (prime Paley orders q0 + 1), two profiles
   construct : round5/construct.md, b = 8,  r = 9M - 8
   lower     : round5/lower.md / round6/sharp.md section 4, b = 33, r = 33(M-1)
and many independent prime choices (consecutive primes = 1 mod 4 from a random start in [1e7, 1e9],
random conjugate per prime), compute EXACT collision moduli m_xy (unit 1, windings k = 0) and the
any-unit moduli, for X = 2M and X = 4M, with D_xy = 2 log m_xy.  Also a variant with iid random primes.

Output: one JSON line per trial in genuine_trials.jsonl; D matrices of the first trials for M <= 32 in
data/.  The aggregation is in genuine_table.py.
Usage: python3 genuine_stats.py [ntrials_construct] [ntrials_lower] [ntrials_iid]
"""
import sys, os, json, math, random, time
import numpy as np
from multiprocessing import Pool
from gauss7 import *
from demands7 import *

MS = [12, 20, 24, 32, 44, 60, 72, 84]
os.makedirs('data', exist_ok=True)
_PROF = {}


def profile(M, variant):
    key = (M, variant)
    if key not in _PROF:
        q0 = M - 1
        S = profile_construct(q0, 8) if variant == 'construct' else profile_lower(q0, 33)
        check_profile(S)
        _PROF[key] = ((1 + S) // 2).astype(np.int64)
    return _PROF[key]


def choose_primes(r, rng, mode):
    if mode == 'window':
        start = rng.randrange(10 ** 7, 10 ** 9)
        return primes_1mod4_window(start, r), start
    ps = set()
    while len(ps) < r:
        n = 4 * rng.randrange(25 * 10 ** 6, 25 * 10 ** 8) + 1      # [1e8, 1e10]
        if is_prime(n):
            ps.add(n)
    return sorted(ps), None


def one_trial(args):
    M, variant, mode, t = args
    Ap = profile(M, variant)
    r = Ap.shape[1]
    rng = random.Random(M * 1000003 + {'construct': 1, 'lower': 2}[variant] * 10007
                        + {'window': 1, 'iid': 2}[mode] * 101 + t)
    ps, start = choose_primes(r, rng, mode)
    gA, gB = gaussian_primes(ps, rng)
    classes = point_classes(Ap, gA, gB, ps, 4 * M)
    rec = dict(M=M, variant=variant, mode=mode, trial=t, r=r, start=start)
    for X in (2 * M, 4 * M):
        m1, many = collision_moduli(classes, X, M)
        for name, mod in (('u1', m1), ('any', many)):
            D = demand_matrix(mod)
            s = stats(D)
            iu = np.triu_indices(M, 1)
            mods = [int(v) for v in mod[iu]]
            s['mmax'] = max(mods)
            # tail of D (pairs with D >= 2 log T)
            s['hist'] = {}
            for v in mods:
                s['hist'][str(v)] = s['hist'].get(str(v), 0) + 1
            rec['%s_X%d' % (name, X // M)] = s
            if M <= 32 and t < 12 and mode == 'window':
                np.save('data/D_%s_M%d_t%d_%s_X%d.npy' % (variant, M, t, name, X // M), D)
                np.save('data/mod_%s_M%d_t%d_%s_X%d.npy' % (variant, M, t, name, X // M),
                        np.array([[int(v) for v in row] for row in mod], dtype=np.int64))
        # window collision counts, primes in (sqrt M, X], level 1, unit 1
        C = window_counts(classes, M, math.sqrt(M), X)
        iu = np.triu_indices(M, 1)
        cnt = np.bincount(C[iu])
        rec['cnt_X%d' % (X // M)] = [int(v) for v in cnt]
    if M <= 32 and t < 12 and mode == 'window':
        with open('data/primes_%s_M%d_t%d.json' % (variant, M, t), 'w') as f:
            json.dump(dict(ps=ps, A=gA, B=gB), f)
    return rec


if __name__ == '__main__':
    nc = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    nl = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    ni = int(sys.argv[3]) if len(sys.argv) > 3 else 60
    jobs = []
    for M in MS:
        jobs += [(M, 'construct', 'window', t) for t in range(nc)]
        jobs += [(M, 'lower', 'window', t) for t in range(nl)]
        jobs += [(M, 'construct', 'iid', t) for t in range(ni)]
    t0 = time.time()
    with Pool(4) as pool, open('genuine_trials.jsonl', 'w') as f:
        for k, rec in enumerate(pool.imap_unordered(one_trial, jobs, chunksize=2)):
            f.write(json.dumps(rec) + '\n')
            if k % 200 == 0:
                print('%d/%d done, %.0fs' % (k, len(jobs), time.time() - t0)); sys.stdout.flush()
    print('done %.0fs' % (time.time() - t0))
