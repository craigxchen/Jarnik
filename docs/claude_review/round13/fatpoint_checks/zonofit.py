# For each certificate printed by prodsearch (lists of balanced z in Z^n), test whether the zonotope
# k0 + sum_i [0,1] z_i' (z' = z without base coordinate) has a translate inside A (i.e. is a genuine zonotope certificate).
import sys, ast, itertools, numpy as np
from fplib import *
fn = sys.argv[1]
cur = None; res = {}
for line in open(fn):
    if line.startswith('=== '):
        name, L, m = line.split()[1:4]; cur = (name, int(L), int(m)); res[cur] = []
    elif line.strip().startswith('[(') and cur:
        res[cur].append(ast.literal_eval(line.strip()))
for (name, L, m), certs in res.items():
    if not certs: continue
    M, pats = MODELS[name](); A = set(box(M, pats, L)); d = M - 1
    nz = 0
    for Z in certs:
        G = np.array([z[1:] for z in Z])
        verts = [G[list(I)].sum(axis=0) if I else np.zeros(d, dtype=int) for r in range(len(G)+1) for I in itertools.combinations(range(len(G)), r)]
        V = np.array(verts); lo = V.min(axis=0); hi = V.max(axis=0)
        ok = any(all(tuple(v + np.array(k0)) in A for v in V) for k0 in itertools.product(*[range(-1-int(l), 2-int(h)) for l, h in zip(lo, hi)]))
        nz += ok
    print(f'{name} L={L} m={m}: certificates listed {len(certs)}, of which zonotopes (translate in A): {nz}')
