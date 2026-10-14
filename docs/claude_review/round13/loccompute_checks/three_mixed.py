"""For every three-cut block (d<=6) where J_w fails: b0 for families from columns u coordinates, and the dimension of the
coordinate support needed."""
import json, re, collections, time
from mixed import *
regions = [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1)]
out = open('three_mixed_out.txt', 'w'); stats = collections.Counter(); t0 = time.time()
for line in open('cuts_three_d6_out.txt'):
    r = json.loads(line[line.index('{'):])
    if r.get('b0') == 1: continue
    sz = eval(re.search(r'regions=(\([^)]*\))', line).group(1)); Lv = eval(re.search(r'L=(\([^)]*\))', line).group(1))
    sets = [[], [], []]; x = 0
    for rg, z in zip(regions, sz):
        for _ in range(z):
            for j in range(3):
                if rg[j]: sets[j].append(x)
            x += 1
    d = x
    res = mixed_b0(*models.cuts(d, sets, Lv), m=r['rho'], bmax=d, verbose=False)
    stats[(r.get('b0'), res['b0_mixed'], d)] += 1
    rec = dict(regions=sz, L=Lv, d=d, rho=r['rho'], firstfail=r['b1_firstfail'], b0_col=r.get('b0'), b0_mixed=res['b0_mixed'])
    out.write(json.dumps(rec) + "\n"); out.flush()
print("(b0_col, b0_mixed, d) -> count:", sorted(stats.items(), key=str), f"[{time.time()-t0:.0f}s]")
open('three_mixed_summary.txt', 'w').write(str(sorted(stats.items(), key=str)) + "\n")
