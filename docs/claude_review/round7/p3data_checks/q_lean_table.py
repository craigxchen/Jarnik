"""Tables of p3data.md §4.2 from q_lean_M{12,20,24}.jsonl, plus the ratio summary including M = 8
(q_results_M8.jsonl, full pipeline).  Usage: python3 q_lean_table.py"""
import json, math, os
import numpy as np

rows = []
for M in (12, 20, 24):
    fn = 'q_lean_M%d.jsonl' % M
    if os.path.exists(fn):
        rows += [json.loads(l) for l in open(fn) if l.strip()]
rows.sort(key=lambda r: (r['M'], r['X'], r['t']))

print('| M | X | t | D_max | 2 lam | W_cut (cert. lb) | W_sc ub (cuts) | W_prime (r, delta) | pairs exact | char4 min | W_cut+sep (greedy) | W_unif |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|')
for r in rows:
    pr = '%.1f (%d, %.1f)' % (r['prime'], r['prime_r'], r['realise_delta']) if r.get('prime') is not None else 'none'
    print('| %d | %d | %d | %.1f | %.1f | %.2f (%.2f) | %.1f (%d) | %s | %s | %s | %.1f | %.1f |' % (
        r['M'], r['X'], r['t'], r['Dmax'], r['lam2'], r['Wcut'], r['Wcut_lb'], r['sc'], r['sc_ncuts'], pr,
        r.get('prime_pairs_ok'), ('%.2f' % r['char4_min']) if r.get('char4_min') is not None else '-',
        r['Wcut_sep_greedy'], r['unif']))

# ratio summary (medians)
m8 = [json.loads(l) for l in open('q_results_M8.jsonl') if l.strip()]
print()
print('| M | X | n | W_cut/2lam | W_sc(ub)/W_cut | W_prime/W_cut | W_prime/2lam | r_prime/M | (W_cut+sep)-W_cut | W_unif/2lam | W_unif/W_prime |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
def med(v):
    v = [x for x in v if x is not None]
    return float(np.median(v)) if v else float('nan')
groups = {}
for r in m8:
    groups.setdefault((8, r['X'] // 8), []).append(dict(lam2=r['lam2'], Wcut=r['Wcut'], sc=r['sc'], prime=r.get('prime'),
                                                    r=r.get('prime_r'), sep=r['Wcut_sep'], unif=r['unif'], M=8))
for r in rows:
    groups.setdefault((r['M'], r['X'] // r['M']), []).append(dict(lam2=r['lam2'], Wcut=r['Wcut'], sc=r['sc'], prime=r.get('prime'),
                                                              r=r.get('prime_r'), sep=r['Wcut_sep_greedy'], unif=r['unif'], M=r['M']))
for (M, xf), g in sorted(groups.items()):
    print('| %d | %dM | %d | %.3f | %.2f | %.2f | %.2f | %.2f | %.1f | %.2f | %.2f |' % (
        M, xf, len(g), med([a['Wcut'] / a['lam2'] for a in g]), med([a['sc'] / a['Wcut'] for a in g]),
        med([a['prime'] / a['Wcut'] if a['prime'] else None for a in g]),
        med([a['prime'] / a['lam2'] if a['prime'] else None for a in g]),
        med([a['r'] / M if a['r'] else None for a in g]), med([a['sep'] - a['Wcut'] for a in g]),
        med([a['unif'] / a['lam2'] for a in g]), med([a['unif'] / a['prime'] if a['prime'] else None for a in g])))
print()
print('note: M = 8 sep column uses the full pipeline positions (placement or equal), M >= 12 the greedy equal spacing')
