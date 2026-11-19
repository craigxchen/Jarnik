"""Structured stress tests: pencils (near-tight BMY), nested trees, long tangency chains."""
import sympy as sp, random
from fractions import Fraction as Fr
from bmy_global import analyse, t
def run(name, cols, N, nl):
    r = analyse(cols, N, nl)
    print(name, {k: (str(v) if isinstance(v, Fr) else v) for k, v in r.items()})
    assert r['slack'] >= 0 and r['iv2'] and r['iv3']
    return r
# pencil y = 1 + u_i t^r : expect slack = 4r
for m in (4, 6, 9):
    for r in (1, 3, 5):
        cols = [(1 + (i+1)*t**r, sp.Integer(1)) for i in range(m)]
        res = run('pencil m=%d r=%d' % (m, r), cols, t**r, [t**r]*m)
        assert res['slack'] == 4*r
# pairs with long tangency at t=0 and t=1, otherwise generic
cols = []
for p in range(4):
    base = p + 1
    for s in (0, 1):
        cols.append((base + (2*s+1)*t**4*(t-1)**3 + p*t, 1 + p*t**2))
run('pairs-chains', cols, (t**4*(t-1)**3), [t**4*(t-1)**3]*8)
# nested tree at 0: groups {0..5} agree mod t, {0,1,2} mod t^3, {0,1} mod t^6
fs = [1 + t + t**3 + t**6*7, 1 + t + t**3 + t**6*11, 1 + t + t**3*5 + t**4, 1 + t*3 + t**2, 1 + t*4 + t**5, 1 + t*9]
fs = [sp.expand(f) for f in fs] + [sp.Integer(2) + t, sp.Integer(5) + t**2]
cols = [(f, sp.Integer(1)) for f in fs]
run('nested', cols, t**6, [t**6, t**6, t**3, t, t, t, sp.Integer(1), sp.Integer(1)])
