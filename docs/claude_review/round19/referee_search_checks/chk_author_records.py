"""Recompute every author record (gain > -40 and all genuinely special) with the referee's own literal formula,
and check the forgetful-pullback data claim (support <= 3 labels, or 4 labels and A/C invariant under x -> 1/x)."""
import sys, json
from collections import Counter
sys.path.insert(0, '.')
from rcore import *
import sympy
AU = '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round19/search_checks/'
xs = sympy.symbols('x0:7')
def forgetful(A, C):
    sup = sorted(set(i for e in A + C for i in e))
    if len(sup) <= 3:
        return 'support%d' % len(sup)
    if len(sup) == 4:
        f = sympy.Mul(*[xs[i] - xs[j] for i, j in A]) / sympy.Mul(*[xs[i] - xs[j] for i, j in C])
        g = f.subs({xs[i]: 1 / xs[i] for i in sup}, simultaneous=True)
        return 'support4-inv' if sympy.simplify(f - g) == 0 else 'support4-NOT'
    return 'support%d' % len(sup)
for D in (4, 5, 6):
    d = json.load(open(AU + 'out/D%d_analysis.json' % D))
    res = Counter(); bad = []
    for r in d['records_gain_gt_-40_and_all_special']:
        A = [tuple(e) for e in r['A']]; C = [tuple(e) for e in r['C']]
        p, q = r['lam']
        g = full_gain(A, C, p, q)
        if g != r['gain']:
            bad.append((r, g))
        res[(g, forgetful(A, C), 'pm1' if (p, q) in ((1, 1), (-1, 1)) else 'special')] += 1
    print(D, 'records', len(d['records_gain_gt_-40_and_all_special']), 'mismatches', len(bad), dict(res))
