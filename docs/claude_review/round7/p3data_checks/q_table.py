"""Tables for p3data.md section 4 from q_results_M*.jsonl."""
import json, glob, math, sys
import numpy as np

rows = []
for f in sorted(glob.glob('q_results_M*.jsonl'), key=lambda s: int(s.split('_M')[1].split('.')[0])):
    rows += [json.loads(l) for l in open(f)]
rows.sort(key=lambda r: (r['M'], r['X'], r['t']))
for r in rows:
    M = r['M']
    D = np.load('data/D_construct_M%d_t%d_u1_X%d.npy' % (M, r['t'], r['X'] // M))
    SM = sum((M - k) * math.log((M - k) / k) for k in range(1, M) if M - k > k)
    r['sep_lb'] = 4.0 / M * (float(D[np.triu_indices(M, 1)].sum()) + 2 * SM)
print('| M | X | t | D_max | 2 lam | W_cut (cert. lb) | W_sc [lb, ub] (cuts) | W_prime (r) | sep LB (Prop 3.2) | W_cut+sep | W_prime+sep (r) | W_unif (= +sep) | char min |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for r in rows:
    sc = '[%.1f, %.1f]%s (%d)' % (r['sc_lb'], r['sc'], '*' if r.get('sc_exact') else '', r['sc_ncuts'])
    pr = '%.1f (%d)' % (r['prime'], r['prime_r']) if r.get('prime') else 'none'
    ps = '%.1f (%d)' % (r['prime_sep'], r['prime_sep_r']) if r.get('prime_sep') else 'none'
    uni = '%.1f%s' % (r['unif'], '' if abs(r['unif_sep'] - r['unif']) < 1e-6 else ' (+sep %.1f)' % r['unif_sep'])
    ok = r.get('prime_pairs_ok', True) and r.get('prime_sep_pairs_ok', True) and r.get('prime_sep_placement_ok', True) and r.get('char_small_ok', True)
    print('| %d | %d | %d | %.1f | %.1f | %.2f (%.2f) | %s | %s | %.1f | %.1f | %s | %s | %s |' % (
        r['M'], r['X'], r['t'], r['Dmax'], r['lam2'], r['Wcut'], r['Wcut_lb'], sc, pr, r.get('sep_lb', float('nan')),
        r['Wcut_sep'], ps, uni, ('%.2f' % r['char_min_all']) if r.get('char_min_all') is not None else '-') + ('' if ok else ' CHECK FAILED'))
print()
print('Summary (medians over instances): ratios')
print('| M | X | n | W_cut/2lam | W_sc(ub)/W_cut | W_prime/W_cut | (W_prime+sep)/(W_cut+sep) | (W_prime+sep)/2lam | W_unif/2lam | W_unif/(W_prime+sep) | r_sep/M |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
from collections import defaultdict
g = defaultdict(list)
for r in rows:
    g[(r['M'], r['X'] // r['M'])].append(r)
for (M, Xf), rs in sorted(g.items()):
    med = lambda f: float(np.median([f(r) for r in rs if f(r) is not None]))
    def safe(f):
        def h(r):
            try:
                return f(r)
            except (KeyError, TypeError):
                return None
        return h
    print('| %d | %dM | %d | %.3f | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f |' % (
        M, Xf, len(rs), med(lambda r: r['Wcut'] / r['lam2']), med(lambda r: r['sc'] / r['Wcut']),
        med(safe(lambda r: r['prime'] / r['Wcut'])), med(safe(lambda r: r['prime_sep'] / r['Wcut_sep'])),
        med(safe(lambda r: r['prime_sep'] / r['lam2'])), med(lambda r: r['unif'] / r['lam2']),
        med(safe(lambda r: r['unif_sep'] / r['prime_sep'])), med(safe(lambda r: r['prime_sep_r'] / M))))
