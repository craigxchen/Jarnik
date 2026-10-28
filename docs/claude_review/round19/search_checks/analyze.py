# Final analysis of a run: full histograms, max gain, exact fibre decompositions of every
# record with gain > -40 (lam = +-1 and special lam), and the forgetful-support test for the gain -32 records.
import json, collections, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-Jarnik/ba4929a4-8665-5790-985a-3f09130a14b9/scratchpad/round19/search_checks')
import exact
fn, outfn = sys.argv[1], sys.argv[2]
hp = collections.Counter(); hg = collections.Counter(); hs = collections.Counter(); st = collections.Counter()
top = []; sp = []; keys = set(); n = 0; tm = 0
for l in open(fn):
    r = json.loads(l); k = (r['orbit'], r['c0'])
    if k in keys: continue
    keys.add(k); n += r['n']; tm += r['time']
    hp.update({int(a): b for a, b in r['hist_pm'].items()}); hg.update({int(a): b for a, b in r['hist_generic'].items()})
    hs.update({int(a): b for a, b in r['hist_special'].items()}); st.update(r['stats']); top += r['top']; sp += r['special']
out = dict(file=fn, chunks=len(keys), pairs=n, cpu_time=round(tm, 1), stats=dict(st),
           hist_lam_pm1=sorted(hp.items(), reverse=True), hist_generic=sorted(hg.items(), reverse=True),
           hist_special_nonsqfree=sorted(hs.items(), reverse=True),
           max_gain=max(list(hp) + list(hg) + list(hs)))
assert not any(isinstance(t[0], str) for t in top), 'mismatch / disc-zero records present'
def lamof(t): return (t[3], 1) if isinstance(t[3], int) else tuple(t[3])
import sympy as sy
xs = sy.symbols('x0:7')
def forgetful(A, C, sup):
    # f = A/C is pulled back from a forgetful map Mbar_{0,8} -> Mbar_{0,4} iff it depends on 4 of the 8 points:
    # support of 3 labels (the 4th point is infinity), or 4 labels and f Moebius-invariant (test: x -> 1/x).
    if len(sup) == 3: return True
    if len(sup) != 4: return False
    fa = sy.Mul(*[xs[i] - xs[j] for i, j in A]); fc = sy.Mul(*[xs[i] - xs[j] for i, j in C])
    inv = {xs[i]: 1 / xs[i] for i in sup}
    return sy.simplify(fa.subs(inv) * fc - fc.subs(inv) * fa) == 0
recs = []
for t in top + [s for s in sp if s[0] <= -40]:
    g, A, C = t[0], [tuple(e) for e in t[1]], [tuple(e) for e in t[2]]
    p, q = lamof(t)
    d = exact.exact_gain(A, C, p, q, detail=True)
    assert d['gain'] == g, (t, d['gain'])
    sup = sorted(set(v for e in A + C for v in e))
    recs.append(dict(gain=g, A=A, C=C, lam=[p, q], support=sup, forgetful=forgetful(A, C, sup), detail=d))
out['records_gain_gt_-40_and_all_special'] = recs
out['n_records'] = len(recs)
out['all_gain_-32_records_are_forgetful_pullbacks'] = all(r['forgetful'] for r in recs if r['gain'] == -32)
out['n_special_nonsqfree'] = len(sp)
out['special_gain_values'] = sorted(collections.Counter(s[0] for s in sp).items(), reverse=True)
json.dump(out, open(outfn, 'w'), indent=1, default=str)
print({k: v for k, v in out.items() if k != 'records_gain_gt_-40_and_all_special' and not k.startswith('hist')})
print('lam=+-1 top', out['hist_lam_pm1'][:6]); print('generic top', out['hist_generic'][:6]); print('special', out['hist_special_nonsqfree'])
