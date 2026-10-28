# Line-independence check: the genuinely special (A, C, lambda) and all histograms must not depend on the lines.
import json, sys, collections
def load(fn):
    sp = set(); hp = collections.Counter(); hg = collections.Counter(); st = collections.Counter()
    for l in open(fn):
        r = json.loads(l)
        for g, A, C, lam in r['special']: sp.add((g, str(A), str(C), tuple(lam)))
        hp.update(r['hist_pm']); hg.update(r['hist_generic']); st.update(r['stats'])
    return sp, hp, hg, st
a = load(sys.argv[1]); b = load(sys.argv[2])
print('special sets identical:', a[0] == b[0], len(a[0]), len(b[0]))
print('lam=+-1 histograms identical:', a[1] == b[1], ' generic identical:', a[2] == b[2])
print('stats A', dict(a[3])); print('stats B', dict(b[3]))
