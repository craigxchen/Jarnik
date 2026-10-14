"""Stronger rigorous LP bound for zeta(A), T = 1, n = 8 models: all ternary functionals f (up to sign) with
b_f = max_A |f.k|, symmetrised over the coordinate permutations that preserve A (verified), dual certificate
rationalised and checked exactly on every candidate v in {-2..2}^d with |f.v| <= 2 b_f for all f.
"""
import sys, itertools, time
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog
from flib import *
from model import profile

name = sys.argv[1]; L = int(sys.argv[2])
M, pats = profile(name)
d = M - 1
A = box(M, pats, L)
Aset = set(A)
Aa = np.array(A, dtype=np.int64)
# automorphisms: coordinate permutations x -> a x + b on F_q (coordinates = F_q, base = inf), plus k -> -k
q = d
perms = []
if name.startswith(('sdc', 'pal')):
    for a in range(1, q):
        for b in range(q):
            perm = [(a * x + b) % q for x in range(q)]
            if all(tuple(k[perm.index(y)] for y in range(q)) in Aset for k in A):
                perms.append(perm)
else:
    perms = [list(range(d))]
print('%s L=%d |A|=%d automorphism perms found: %d' % (name, L, len(A), len(perms)), flush=True)
F = [f for f in itertools.product((-1, 0, 1), repeat=d) if any(f) and next(x for x in f if x) == 1]
Fa = np.array(F, dtype=np.int64)
bf = np.abs(Aa @ Fa.T).max(axis=0)
V = np.array([v for v in itertools.product(range(-2, 3), repeat=d) if any(v) and next(x for x in v if x) > 0], dtype=np.int64)
FV = np.abs(V @ Fa.T)
ok = (FV <= 2 * bf[None, :]).all(axis=1)
V = V[ok]; FV = FV[ok]
print('functionals %d, candidate directions %d' % (len(F), len(V)), flush=True)


def canon_f(f):
    f = tuple(f)
    if next(x for x in f if x) < 0:
        f = tuple(-x for x in f)
    return f


fidx = {f: i for i, f in enumerate(F)}
# orbits of functionals under perms (and sign)
orb_id = [-1] * len(F)
orbs = []
for i, f in enumerate(F):
    if orb_id[i] >= 0:
        continue
    o = set()
    for p in perms:
        g = [0] * d
        for x in range(d):
            g[p[x]] = f[x]
        o.add(fidx[canon_f(g)])
    for j in o:
        orb_id[j] = len(orbs)
    orbs.append(sorted(o))
# constraint matrix: per v, per functional-orbit, sum of |f.v|
G = np.zeros((len(V), len(orbs)))
for oi, o in enumerate(orbs):
    G[:, oi] = FV[:, o].sum(axis=1)
G, inv = np.unique(G, axis=0, return_inverse=True)
caps = np.array([2 * bf[o].sum() for o in orbs], dtype=float)
res = linprog(caps, A_ub=-G, b_ub=-np.ones(len(G)), bounds=[(0, None)] * len(orbs), method='highs')
print('LP value %.6f (%d orbit variables, %d distinct constraints)' % (res.fun, len(orbs), len(G)), flush=True)
y = [Fraction(x).limit_denominator(10 ** 5) for x in res.x]
Gi = G.astype(np.int64)
worst = min(sum(Fraction(int(a)) * b for a, b in zip(row, y)) for row in Gi)
if worst < 1:
    y = [b / worst for b in y]
assert all(sum(Fraction(int(a)) * b for a, b in zip(row, y)) >= 1 for row in Gi)
bound = sum(Fraction(int(2 * bf[o].sum())) * b for o, b in zip(orbs, y))
print('RESULT %s L=%d: rigorous zeta <= %s = %.5f' % (name, L, bound, float(bound)))
