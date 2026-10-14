"""b0 versus ratio rho/(2 L_max) for three-cut blocks (d<=6) and two-cut blocks (d<=7)."""
import json, re, collections
def ratio_bin(x): return f"{int(x*2)/2:.1f}-{int(x*2)/2+0.5:.1f}"
for fn, kind in (('cuts_two_d7_out.txt', 'two'), ('cuts_three_d6_out.txt', 'three')):
    mixed = {}
    if kind == 'three':
        for l in open('three_mixed_out.txt'):
            r = json.loads(l); mixed[(tuple(r['regions']), tuple(r['L']))] = r['b0_mixed']
    tab = collections.defaultdict(collections.Counter)
    for line in open(fn):
        r = json.loads(line[line.index('{'):])
        Lv = eval(re.search(r'L=(\([^)]*\))', line).group(1))
        rat = r['rho'] / (2 * max(Lv))
        key = r.get('b0')
        if kind == 'three':
            sz = eval(re.search(r'regions=(\([^)]*\))', line).group(1))
            if key != 1: key = f"col={key},mixed={mixed.get((tuple(sz), tuple(Lv)), '?')}"
        tab[ratio_bin(rat)][key] += 1
    print(kind)
    for k in sorted(tab): print("  rho/(2Lmax) in", k, dict(tab[k]))
