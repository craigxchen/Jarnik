"""Character margins of the realised actual-prime profiles of §4.2 (ils2_M*.jsonl, best_ub_M*.jsonl).
For each profile: re-check the pair inequalities exactly (verify_profile.check_pairs), the exact
support-<= 4 character check (verify_profile.char_check_small: integers) when M <= 24, and
min s_c over all zero-sum integer c != 0 with |c_x| <= 12 (verify_profile.char_min_milp: M MILPs,
floating point), reporting whether the minimiser is a pair.  Output: char_all_out.txt / .json."""
import sys, os, json, glob, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); UP = os.path.dirname(HERE)
sys.path.insert(0, UP)
from verify_profile import check_pairs, char_check_small, char_min_milp   # noqa: E402

pats = sys.argv[1:] or ['ils2_M*.jsonl', 'best_ub_M*.jsonl']
KEY = {'ils2': 'prime_ils_cols', 'best': 'prime_cols'}
out = []
for fn in pats:
    key = KEY[fn.split('_')[0]]
    for f in sorted(glob.glob(os.path.join(HERE, fn))):
        for line in open(f):
            r = json.loads(line)
            if r['M'] == 8 or key not in r:
                continue
            M, X, t = r['M'], r['X'], r['t']
            cols = [(m, pl) for m, pl in r[key]]
            mod = np.load(os.path.join(UP, 'data/mod_construct_M%d_t%d_u1_X%d.npy' % (M, t, X // M)))
            ok, worst = check_pairs(cols, mod, X)
            t0 = time.time()
            small = char_check_small(cols, M, X) if M <= 24 else (None, None)
            cmin, carg = char_min_milp(cols, M, bound=12, time_limit=60)
            supp = int((np.abs(carg) > 0).sum()) if carg is not None else None
            rec = dict(src=os.path.basename(f), M=M, X=X, t=t, W=float(sum(np.log(p) for _, pl in cols for p in pl)),
                       pairs_ok=bool(ok), small_ok=small[0], small_min=(small[1][1] if small[0] else None),
                       char_min=cmin, argmin_support=supp, time=time.time() - t0)
            out.append(rec)
            print('%s M=%d X=%d t=%d W=%.1f pairs_ok=%s supp4_exact_ok=%s (min %.2f) | min s_c (|c|<=12) = %.3f at support %s  %.0fs' % (
                rec['src'], M, X, t, rec['W'], ok, small[0], rec['small_min'] if rec['small_min'] is not None else float('nan'),
                cmin, supp, rec['time']), flush=True)
tag = '_'.join(x.replace('*', '').replace('.jsonl', '') for x in pats)
json.dump(out, open(os.path.join(HERE, 'char_all_%s.json' % tag), 'w'), indent=0)
