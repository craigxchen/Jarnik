"""Markdown trend tables for compute.md §6 from ../boxes.csv (run table.py first)."""
import csv
from collections import defaultdict
B = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/'
rows = [r for r in csv.DictReader(open(B + 'boxes.csv')) if r.get('zeta_all') not in (None, '')]
def rp(r):
    return int(r['reach_pair_P1'] if r['reach_pair_P1'] not in ('', None) else r['reach_pair_P2'])
def best(r):
    return max(int(r[k]) for k in ('reach_pair_P1', 'reach_cube_P1', 'reach_pair_P2', 'reach_cube_P2') if r.get(k) not in ('', None))
def stats(rs):
    z = [int(r['zeta_all']) for r in rs]
    p = [rp(r) for r in rs]; b = [best(r) for r in rs]
    rho = [int(r['rho']) for r in rs]
    f = [i for i in range(len(rs)) if z[i] < rho[i]]
    return dict(n=len(rs), gt=sum(p[i] > z[i] for i in range(len(rs))), eq=sum(p[i] == z[i] for i in range(len(rs))),
                lt=sum(p[i] < z[i] for i in range(len(rs))), bgt=sum(b[i] > z[i] for i in range(len(rs))),
                beq=sum(b[i] == z[i] for i in range(len(rs))), blt=sum(b[i] < z[i] for i in range(len(rs))),
                mgap=sum(p[i] - z[i] for i in range(len(rs))) / len(rs), mbgap=sum(b[i] - z[i] for i in range(len(rs))) / len(rs),
                mrr=sum(p[i] / rho[i] for i in range(len(rs))) / len(rs), mzr=sum(z[i] / rho[i] for i in range(len(rs))) / len(rs),
                mrl=sum(float(r['rho_over_L']) for r in rs) / len(rs), nf=len(f), fgt=sum(b[i] > z[i] for i in f),
                share=sum(b[i] > z[i] for i in range(len(rs))) / len(rs))
def table(title, groups, order=None):
    print('\n' + title + '\n')
    print('| group | boxes | faithful (`zeta < rho`) | pair: > / = / < `zeta` | best (pair or cube): > / = / < | share best > `zeta` | mean reach_pair − `zeta` | mean best − `zeta` | mean reach_pair/`rho` | mean `zeta`/`rho` | mean `rho/L` |')
    print('| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |')
    for k in (order or sorted(groups)):
        if k not in groups: continue
        s = stats(groups[k])
        print('| %s | %d | %d | %d / %d / %d | %d / %d / %d | %.2f | %+.2f | %+.2f | %.3f | %.3f | %.2f |' % (
            k, s['n'], s['nf'], s['gt'], s['eq'], s['lt'], s['bgt'], s['beq'], s['blt'], s['share'], s['mgap'], s['mbgap'],
            s['mrr'], s['mzr'], s['mrl']))
def G(key):
    g = defaultdict(list)
    for r in rows:
        k = key(r)
        if k is not None: g[k].append(r)
    return g
israndom = lambda r: r['family'].startswith('random')
table('**Table 6.1 (by `d`, all boxes).**', G(lambda r: 'd = %s' % r['d']))
table('**Table 6.2 (by `d` and kind).**', G(lambda r: 'd = %s, %s' % (r['d'], 'random' if israndom(r) else 'structured')))
table('**Table 6.3 (faithful boxes only, `zeta < rho`, by `d`).**', G(lambda r: ('d = %s' % r['d']) if int(r['zeta_all']) < int(r['rho']) else None))
def rl(r):
    x = float(r['rho_over_L'])
    return 'rho/L <= 3' if x <= 3 else ('3 < rho/L <= 4' if x <= 4 else ('4 < rho/L <= 5' if x <= 5 else ('5 < rho/L <= 6' if x <= 6 else 'rho/L > 6')))
table('**Table 6.4 (by `rho/L`, all boxes).**', G(rl), ['rho/L <= 3', '3 < rho/L <= 4', '4 < rho/L <= 5', '5 < rho/L <= 6', 'rho/L > 6'])
def sz(r):
    s = int(r['size'])
    return '|A| < 100' if s < 100 else ('100 <= |A| < 500' if s < 500 else ('500 <= |A| < 2000' if s < 2000 else '|A| >= 2000'))
table('**Table 6.5 (by box size `|A|`, all boxes).**', G(sz), ['|A| < 100', '100 <= |A| < 500', '500 <= |A| < 2000', '|A| >= 2000'])
table('**Table 6.6 (random boxes by `d` and `L`).**', G(lambda r: ('d = %s, L = %s' % (r['d'], r['L'])) if israndom(r) else None))
