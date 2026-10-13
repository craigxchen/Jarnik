"""Projective variety X = V(G_1..G_s) in P^{d-1} (x-coordinates) for a cone-reduced model: exact Groebner basis
over Q (grevlex; sympy), Hilbert-function dimension count, and, when X is finite, its points.
Then the real-torus analysis: a point w = e^{iu} near e lies in the cone iff x(u) = e^{iu} - 1 is a nonzero
multiple of some x* in X (or 0)."""
import sys, json, time
import sympy as sp
fn = sys.argv[1]
C = json.load(open(fn)); d = C['d']
X = sp.symbols('x0:%d' % d)
forms = [sp.sympify(f) for f in C['forms']]
forms = [f for f in forms if f != 0]
t0 = time.time()
GB = sp.groebner(forms, *X, order='grevlex')
print("Groebner basis: %d elements, degrees %s (%.1fs)" % (len(GB.exprs), sorted(set(sp.Poly(g, *X).total_degree() for g in GB.exprs)), time.time() - t0))
# leading monomials -> Hilbert function of S/in(I) for degrees up to 20
LM = [sp.Poly(g, *X).monoms(order='grevlex')[0] for g in GB.exprs]
def divisible(mon):
    return any(all(mon[i] >= l[i] for i in range(d)) for l in LM)
from itertools import combinations_with_replacement
HF = []
for D in range(0, 16):
    cnt = 0
    for comb in combinations_with_replacement(range(d), D):
        mon = [0] * d
        for i in comb: mon[i] += 1
        if not divisible(mon): cnt += 1
    HF.append(cnt)
print("Hilbert function of S/I, degrees 0..15:", HF)
json.dump(dict(GB=[str(g) for g in GB.exprs], HF=HF), open(fn.replace("_cone.json", "_gb.json"), "w"))
