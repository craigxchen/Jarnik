"""Tables 4.2a and 4.2b of p3data.md from: ils2_M*.jsonl, best_ub_M*.jsonl, ../q_lean_M*.jsonl,
eqlb_out.json, wcut_extra.json, ../perron_check.json, ../q_results_M8.jsonl (M = 8, §4.1).
Best W_sc upper bound = min over methods; best W_prime = min over verified actual-prime profiles.
Usage: python3 table42.py > table42_out.txt"""
import json, os, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); UP = os.path.dirname(HERE)


def jl(path):
    out = []
    for f in glob.glob(path):
        out += [json.loads(l) for l in open(f) if l.strip()]
    return out


eq = {(d['M'], d['X'], d['t']): d for d in json.load(open(os.path.join(HERE, 'eqlb_out.json')))}
lean = {(r['M'], r['X'], r['t']): r for r in jl(os.path.join(UP, 'q_lean_M*.jsonl'))}
ils = {(r['M'], r['X'], r['t']): r for r in jl(os.path.join(HERE, 'ils2_M*.jsonl')) if r['M'] > 8}
bub = {(r['M'], r['X'], r['t']): r for r in jl(os.path.join(HERE, 'best_ub_M*.jsonl')) if r['M'] > 8}
cert = {(r['M'], r['X'], r['t']): (r['W'], r['lb']) for r in json.load(open(os.path.join(UP, 'perron_check.json')))}
if os.path.exists(os.path.join(HERE, 'wcut_extra.json')):
    for r in json.load(open(os.path.join(HERE, 'wcut_extra.json'))):
        cert[(r['M'], r['X'], r['t'])] = (r['W'], r['lb'])
for k, r in lean.items():
    cert.setdefault(k, (r['Wcut'], r['Wcut_lb']))
m8 = {(r['M'], r['X'], r['t']): r for r in jl(os.path.join(UP, 'q_results_M8.jsonl'))}

keys = sorted(set(ils) | set(bub))
rows = []
for k in keys:
    M, X, t = k
    a, b, ln = ils.get(k, {}), bub.get(k, {}), lean.get(k, {})
    base = a or b
    wc = cert.get(k, (None, None))
    ubs = {'ils': a.get('W_ils'), 'sa': b.get('W_ub'), 'rwl1': ln.get('sc')}
    ubs = {m: v for m, v in ubs.items() if v is not None}
    best_m = min(ubs, key=ubs.get)
    pr = {'ils': (a.get('prime_ils'), a.get('prime_ils_r'), a.get('stair'), a.get('char4_min')),
          'sa': (b.get('prime'), b.get('prime_r'), b.get('stair'), b.get('char4_min')),
          'rwl1': (ln.get('prime'), ln.get('prime_r'), ln.get('stair'), ln.get('char4_min'))}
    pr = {m: v for m, v in pr.items() if v[0] is not None}
    bp = min(pr, key=lambda m: pr[m][0]) if pr else None
    r_ub = {'ils': a.get('ils_r'), 'sa': b.get('ub_r'), 'rwl1': ln.get('sc_ncuts')}[best_m]
    rows.append(dict(M=M, X=X, t=t, Dmax=base['Dmax'], lam2=base['lam2'], Wcut=wc[0], Wlb=wc[1],
                     LB_eq=eq[k]['LB_eq'], ubs=ubs, best=ubs[best_m], best_m=best_m, best_r=r_ub,
                     prime=pr[bp][0] if bp else None, prime_r=pr[bp][1] if bp else None,
                     stair=pr[bp][2] if bp else None, char4=pr[bp][3] if bp else None, prime_m=bp,
                     unif=2 * (M - 1) * base['Dmax']))
for k, r in sorted(m8.items()):
    rows.append(dict(M=8, X=k[1], t=k[2], Dmax=r['Dmax'], lam2=r['lam2'], Wcut=r['Wcut'], Wlb=r['Wcut_lb'],
                     LB_eq=eq[k]['LB_eq'], ubs={'milp': r['sc']}, best=r['sc'], best_m='milp*' if r.get('sc_exact') else 'milp',
                     best_r=r['sc_ncuts'], prime=r['prime'], prime_r=r['prime_r'], stair=None,
                     char4=r.get('char_small_min'), prime_m='q_run', unif=r['unif']))
rows.sort(key=lambda r: (r['M'], r['X'], r['t']))

gaps = [abs(r['Wcut'] - r['Wlb']) / r['Wcut'] for r in rows if r['Wcut'] and r['Wlb']]
print('largest relative gap W_cut vs certified dual bound over all rows: %.1e' % max(gaps))
print()
f = lambda v: '--' if v is None else '%.1f' % v
print('Table 4.2a (M >= 12)')
print('| M | X | t | D_max | 2 lam | W_cut | LB_eq | W_ils | W_sa | rwl1 | best (method, r) | best/2lam | W_prime (r) | stair(r) | W_unif | char4 min |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
for r in rows:
    if r['M'] == 8:
        continue
    u = r['ubs']
    print('| %d | %d | %d | %.1f | %.1f | %s | %s | %s | %s | %s | %.1f (%s, %s) | %.2f | %s (%s) | %s | %.1f | %s |' % (
        r['M'], r['X'], r['t'], r['Dmax'], r['lam2'], f(r['Wcut']), f(r['LB_eq']), f(u.get('ils')), f(u.get('sa')),
        f(u.get('rwl1')), r['best'], r['best_m'], r['best_r'], r['best'] / r['lam2'], f(r['prime']), r['prime_r'],
        f(r['stair']), r['unif'], f(r['char4'])))
print()
print('Table 4.2b: ratios to 2 lambda_max, median [min, max]')
print('| M | X | n | W_cut | LB_eq | best W_sc UB | W_prime | W_unif | W_unif / best UB | best r / M |')
print('|---|---|---|---|---|---|---|---|---|---|')
for M in sorted({r['M'] for r in rows}):
    for X in sorted({r['X'] for r in rows if r['M'] == M}):
        g = [r for r in rows if r['M'] == M and r['X'] == X]

        def st(v):
            v = [x for x in v if x is not None]
            return '%.2f [%.2f, %.2f]' % (np.median(v), min(v), max(v)) if v else '--'
        print('| %d | %d | %d | %s | %s | %s | %s | %s | %s | %s |' % (
            M, X, len(g), st([r['Wcut'] / r['lam2'] if r['Wcut'] else None for r in g]),
            st([r['LB_eq'] / r['lam2'] for r in g]), st([r['best'] / r['lam2'] for r in g]),
            st([r['prime'] / r['lam2'] if r['prime'] else None for r in g]),
            st([r['unif'] / r['lam2'] for r in g]), st([r['unif'] / r['best'] for r in g]),
            st([r['best_r'] / M for r in g])))
json.dump(rows, open(os.path.join(HERE, 'table42_rows.json'), 'w'), indent=0)
