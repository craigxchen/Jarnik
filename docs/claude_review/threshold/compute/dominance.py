"""Empirical check: is lam -> reg_lam(A) order-reversing for the dominance order?"""
import json, glob, sys
from modla import P1
kind = sys.argv[1] if len(sys.argv) > 1 else 'shell'
def dominates(a, b):  # a >= b in dominance
    sa = sb = 0
    for i in range(max(len(a), len(b))):
        sa += a[i] if i < len(a) else 0
        sb += b[i] if i < len(b) else 0
        if sa < sb: return False
    return True
tot = viol = 0; ex = []
for fn in sorted(glob.glob(f'data/{kind}_M*.jsonl')):
    for line in open(fn):
        r = json.loads(line)
        if r['prime'] != P1: continue
        C = [(tuple(c['lam']), c['reg']) for c in r['comps']]
        for a, ra in C:
            for b, rb in C:
                if a != b and dominates(a, b):
                    tot += 1
                    if ra > rb:
                        viol += 1
                        if len(ex) < 15: ex.append((r['M'], r['s'], r['t'], a, ra, b, rb))
print('dominance pairs', tot, 'violations (lam >= mu but reg_lam > reg_mu)', viol)
for e in ex: print(e)
