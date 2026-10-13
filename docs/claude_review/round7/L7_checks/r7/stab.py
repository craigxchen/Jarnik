"""EXACT: stabiliser in S_7 of the cub13 D-vector."""
import itertools
from fibcone import vec, norm, planes13
D = vec(3, mabc={p: 1 for p in planes13})
Dn = {}
for k, v in D.items(): Dn[norm(k)] = Dn.get(norm(k), 0) + v
stab = []
for perm in itertools.permutations(range(1, 8)):
    sg = dict(zip(range(1, 8), perm))
    E = {}
    for k, v in Dn.items():
        kk = norm({sg[x] for x in k}); E[kk] = E.get(kk, 0) + v
    if E == Dn: stab.append(perm)
print('stabiliser order', len(stab))
for p in stab: print([ (i, p[i-1]) for i in range(1, 8) if p[i-1] != i])
