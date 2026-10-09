"""Independent exact construction of the generalized alternant on S^5(6,3) (upper.md Remark 5.4) and
its forced-divisibility ratio ord / MW(supp)."""
import itertools
from fractions import Fraction as Fr
from verify_B import MW
reps = [(5,1,0,-1,-2), (4,2,0,-1,-2), (3,2,1,0,-3), (3,2,1,-1,-2)]
def vand(v):
    p = 1
    for i in range(5):
        for j in range(i+1, 5): p *= (v[i] - v[j])
    return p
def esym(v, r):
    return sum(Fr(1) * eval('*'.join(str(v[i]) for i in c)) if c else 1 for c in itertools.combinations(range(5), r))
rows = [[vand(v) * esym(v, r) for v in reps] for r in (0, 2, 3)]
# kernel of 3x4 integer matrix
import numpy as np
A = [[Fr(x) for x in r] for r in rows]
# solve by brute force small integers
best = None
for g in itertools.product(range(-6, 7), repeat=4):
    if all(x == 0 for x in g): continue
    if all(sum(a * b for a, b in zip(r, g)) == 0 for r in A):
        best = g; break
print("gamma =", best)
def sgn(p):
    s = 1
    for i in range(len(p)):
        for j in range(i+1, len(p)):
            if p[i] > p[j]: s = -s
    return s
c = {}
for g, v in zip(best, reps):
    for p in itertools.permutations(range(5)):
        k = tuple(v[p[i]] for i in range(5))
        c[k] = c.get(k, 0) + g * sgn(p)
c = {k: x for k, x in c.items() if x}
print("support size", len(c), " shell (p,n):", set((sum(max(x,0) for x in k), sum(max(-x,0) for x in k)) for k in c))
# order: all moments of degree < m vanish (exact, in the first 4 coordinates)
def order(c, maxd=20):
    pts = list(c.keys())
    for d in range(maxd):
        for e in itertools.product(range(d+1), repeat=4):
            if sum(e) != d: continue
            mom = sum(cv * k[0]**e[0] * k[1]**e[1] * k[2]**e[2] * k[3]**e[3] for k, cv in c.items())
            if mom != 0: return d
    return None
m = order(c)
S = [k[:4] for k in c]
mw = MW(S, 4)
print("ord =", m, " t = p+n =", 9, " ord/t =", Fr(m, 9), " MW(supp) =", mw, " ord/MW =", Fr(m) / mw)
