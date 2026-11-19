"""Take a strict integer solution of the k=1,2 Im-moment conditions (found by search), G_e = T + c_e with T a
large real integer.  Prediction: cofactors are cubic polynomials in T, Delta ~ P^-3 (nu = 3); the division
step (VI.6) still applies to pairs with m = 2 (needs only nu < 2m) and forces Im Q(T) = 0.  Moving B off the
real line by a Gaussian offset destroys this (nu jumps to ~5)."""
import itertools, math, sys
import numpy as np, sympy as sp
from itertools import combinations
sys.path.insert(0, '.')
from gauss import *
edges = EDGES
Inc = np.array([[int(i in e) for e in edges] for i in range(5)]); E = Inc[:4]-Inc[4]
def imstars(x, y, k):
    v = [sp.im(sp.expand((sp.Integer(a)+sp.I*b)**k)) for a, b in zip(x, y)]
    return [sum(v[e] for e in range(10) if i in edges[e]) for i in range(5)]
# small search for one strict example
found = None
rng = [v for v in range(-3, 4) if v]
for yf in itertools.product(rng, repeat=6):
    y = np.zeros(10, dtype=int); y[:6] = yf
    sol = np.linalg.lstsq(E[:, 6:].astype(float), -(E[:, :6] @ np.array(yf)).astype(float), rcond=None)[0]
    if np.abs(sol-np.rint(sol)).max() > 1e-9: continue
    y[6:] = np.rint(sol).astype(int)
    if (y == 0).any() or np.abs(y).max() > 3: continue
    M = sp.Matrix(E*y[None, :]); ns = M.nullspace()
    for co in itertools.product(range(-2, 3), repeat=len(ns)):
        xv = sum((cc*v for cc, v in zip(co, ns)), sp.zeros(10, 1))
        if any(sp.fraction(q)[1] != 1 for q in xv): continue
        x = [int(q) for q in xv]
        if max(map(abs, x)) > 6: continue
        if len(set(zip(x, y.tolist())) | set(zip(x, (-y).tolist()))) == 20:
            found = (x, y.tolist()); break
    if found: break
x, y = found
print('example c_e = x+iy:', list(zip(x, y)))
print('Im star sums k=1..3:', [imstars(x, y, k) for k in (1, 2, 3)])
for T in [10**4, 10**8, 10**16]:
    G = [(T+a, b) for a, b in zip(x, y)]
    z = tuple_from_blocks(G); H, logD = width_and_H(z)
    print(f'T=1e{int(math.log10(T))}: nu = {6+20*logD/math.log(H):.3f}')
    Gs = [(T+a, b+1) for a, b in zip(x, y)]   # B = T + i, off the line
    z = tuple_from_blocks(Gs); H, logD = width_and_H(z)
    print(f'      B=T+i: nu = {6+20*logD/math.log(H):.3f}')
# division check: pairs with m=2, Im Q(T) for real T should vanish identically (Q real)
u = sp.Symbol('u'); roots = [-(a+sp.I*b) for a, b in zip(x, y)]
for (i, j) in combinations(range(5), 2):
    Ai = [k for k, e in enumerate(edges) if i in e and e != (i, j)]
    Dj = [k for k, e in enumerate(edges) if j in e and e != (i, j)]
    A = sp.expand(sp.prod([u-roots[k] for k in Ai])); D = sp.expand(sp.prod([u-roots[k] for k in Dj]))
    diff = sp.Poly(sp.expand(A-D), u)
    if diff.degree() == 2:
        Q, R = sp.div(sp.Poly(D, u), diff)
        T = 10**8
        print('pair', (i, j), 'm=2, Q =', Q.as_expr(), ' Im Q(1e8) =', sp.im(sp.expand(Q.as_expr().subs(u, T))))
