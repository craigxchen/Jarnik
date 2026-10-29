"""Trend tables from ../boxes.csv: share of reach > / = / < zeta, mean reach - zeta, reach/rho, by d, by family,
by rho/L bin and by |A| bin.  Writes ../trend.csv."""
import csv
from collections import defaultdict
B = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round28/'
rows = [r for r in csv.DictReader(open(B + 'boxes.csv')) if r.get('zeta_all') not in (None, '')]
def rp(r):
    return int(r['reach_pair_P1'] if r['reach_pair_P1'] not in ('', None) else r['reach_pair_P2'])
def best_reach(r):
    v = [int(r[k]) for k in ('reach_pair_P1', 'reach_cube_P1', 'reach_pair_P2', 'reach_cube_P2') if r.get(k) not in ('', None)]
    return max(v)
def grp(key):
    g = defaultdict(list)
    for r in rows:
        g[key(r)].append(r)
    return g
def summarize(name, g):
    out = []
    for k in sorted(g, key=lambda x: (str(type(x)), x)):
        rs = g[k]
        dp = [rp(r) - int(r['zeta_all']) for r in rs]
        dc = [best_reach(r) - int(r['zeta_all']) for r in rs]
        rr = [rp(r) / int(r['rho']) for r in rs]
        zr = [int(r['zeta_all']) / int(r['rho']) for r in rs]
        out.append(dict(grouping=name, group=k, boxes=len(rs),
                        pair_gt=sum(x > 0 for x in dp), pair_eq=sum(x == 0 for x in dp), pair_lt=sum(x < 0 for x in dp),
                        best_gt=sum(x > 0 for x in dc), best_eq=sum(x == 0 for x in dc), best_lt=sum(x < 0 for x in dc),
                        mean_reach_minus_zeta=round(sum(dp) / len(dp), 3), mean_best_minus_zeta=round(sum(dc) / len(dc), 3),
                        mean_reach_over_rho=round(sum(rr) / len(rr), 3), mean_zeta_over_rho=round(sum(zr) / len(zr), 3),
                        mean_rho_over_L=round(sum(float(r['rho_over_L']) for r in rs) / len(rs), 3),
                        faithful=sum(int(r['zeta_all']) < int(r['rho']) for r in rs),
                        faithful_pair_gt=sum(int(r['zeta_all']) < int(r['rho']) and rp(r) > int(r['zeta_all']) for r in rs),
                        faithful_best_gt=sum(int(r['zeta_all']) < int(r['rho']) and best_reach(r) > int(r['zeta_all']) for r in rs)))
    return out
def rhoL_bin(r):
    x = float(r['rho_over_L'])
    return '<=3' if x <= 3 else ('(3,4]' if x <= 4 else ('(4,5]' if x <= 5 else ('(5,6]' if x <= 6 else '>6')))
def size_bin(r):
    s = int(r['size'])
    return '<100' if s < 100 else ('100-499' if s < 500 else ('500-1999' if s < 2000 else '>=2000'))
allrows = []
allrows += summarize('d', grp(lambda r: int(r['d'])))
allrows += summarize('family', grp(lambda r: r['family']))
allrows += summarize('d,family', grp(lambda r: '%s,%s' % (r['d'], 'random' if r['family'].startswith('random') else 'structured')))
allrows += summarize('rho/L', grp(rhoL_bin))
allrows += summarize('size', grp(size_bin))
allrows += summarize('L', grp(lambda r: int(r['L'])))
allrows += summarize('faithful', grp(lambda r: 'zeta<rho' if int(r['zeta_all']) < int(r['rho']) else 'zeta=rho'))
allrows += summarize('d,L(random)', grp(lambda r: '%s,L=%s,%s' % (r['d'], r['L'], 'rand' if r['family'].startswith('random') else 'str')))
with open(B + 'trend.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(allrows[0].keys())); w.writeheader()
    for r in allrows: w.writerow(r)
for r in allrows:
    print("%-12s %-22s n=%3d pair >/=/< %2d/%2d/%2d  best >/=/< %2d/%2d/%2d  mean(r-z)=%6.2f  reach/rho=%.3f zeta/rho=%.3f rho/L=%.2f faithful=%d (crit>zeta among them: pair %d, best %d)" % (
        r['grouping'], r['group'], r['boxes'], r['pair_gt'], r['pair_eq'], r['pair_lt'], r['best_gt'], r['best_eq'], r['best_lt'],
        r['mean_reach_minus_zeta'], r['mean_reach_over_rho'], r['mean_zeta_over_rho'], r['mean_rho_over_L'], r['faithful'], r['faithful_pair_gt'], r['faithful_best_gt']))
