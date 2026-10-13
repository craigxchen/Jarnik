import numpy as np
import cub13f as F
from collections import Counter
P = np.polynomial.polynomial
rng = np.random.default_rng(2)
base = dict(tau=-2.0+0j, s34=3.0+0j, s35=2.5+0j, s45=2.457710+0j)
names = ['tau', 's34', 's35', 's45']
def res(v, var):
    a = dict(base); a[var] = v[0]
    d = F.build(a['tau'], a['s34'], a['s35'], a['s45'])
    return np.array([P.polyval(v[1], d['R']), P.polyval(v[1], d['T'])]), d
def newton(v, var):
    v = np.array(v, complex)
    for _ in range(80):
        f, _ = res(v, var); h = 1e-7; J = np.zeros((2, 2), complex)
        for k in range(2):
            dv = np.zeros(2, complex); dv[k] = h
            J[:, k] = (res(v + dv, var)[0] - res(v - dv, var)[0]) / (2 * h)
        try: st = np.linalg.solve(J, -f)
        except Exception: return None
        if not np.all(np.isfinite(st)): return None
        v = v + st
        if np.linalg.norm(st) < 1e-13 * (1 + np.linalg.norm(v)): break
    f, d = res(v, var)
    if not np.all(np.isfinite(f)) or np.linalg.norm(f) > 1e-9: return None
    return v, d
def ok(d):
    cnt = Counter()
    for p, big in F.events(d['c']):
        for g in big: cnt[g] += 1
    trip = {g for g in cnt if len(g) == 3}
    pairs = {g: cnt[g] for g in cnt if len(g) == 2}
    return trip == F.TRIPLES13 and all(cnt[g] == 1 for g in trip) and pairs == {(1,4):1,(1,5):1,(2,3):1,(1,6):3} and not any(len(g) >= 4 for g in cnt)
for var in names:
    if var == 's45':
        base['s45'] = None
    sols = []
    for k in range(600):
        v0 = rng.normal(size=2) * 3 + 1j * rng.normal(size=2) * 3
        r = newton(v0, var) if var != 's45' or True else None
        if r is None: continue
        v, d = r
        if not ok(d): continue
        if all(abs(v[0] - w[0]) > 1e-6 for w in sols): sols.append(v)
    print(var, 'genuine solutions:', len(sols), [np.round(w[0], 5) for w in sols])
    if var == 's45': base['s45'] = 2.457710+0j
