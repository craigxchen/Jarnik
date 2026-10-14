# validation: find lattice zonotopes sum_i [0,1] g_i (+k0) inside A, test membership of prod <g_i,u> in T_m
import sys, itertools, numpy as np
from fplib import *
name, L, m = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
M, pats = MODELS[name](); A = box(M, pats, L); d = M-1
Aset = set(A)
monos_m, Tb = tangent_space_mod(A, m)
idx = {a:i for i,a in enumerate(monos_m)}; rT = Tb.shape[0]
def in_T(Z):
    poly = product_coeffs(Z, m)
    vec = [0]*len(monos_m)
    for e,c in poly.items():
        if max(e) > 2: return False
        vec[idx[e]] = (c % PR) * pow(multinom(m,e) % PR, -1, PR) % PR
    rows = Tb.tolist() + [vec]
    return nmod_mat(len(rows), len(monos_m), [x for r in rows for x in r], PR).rank() == rT
# generators: k'-vectors with entries in {-1,0,1}, small support
gens = [g for g in itertools.product((-1,0,1), repeat=d) if 0 < sum(map(abs,g)) <= 2]
gens = [g for g in gens if next(x for x in g if x) > 0]
rng = np.random.default_rng(0)
found = 0; tested = 0
for trial in range(20000):
    G = [gens[i] for i in rng.integers(len(gens), size=m)]
    # zonotope vertices relative to k0; choose k0 = -sum of generators with negative... try k0 = -(sum G)//2 rounded variants
    S = np.array(G)
    verts = [tuple(S[list(I)].sum(axis=0)) if I else tuple([0]*d) for r in range(m+1) for I in itertools.combinations(range(m), r)]
    V = np.array(verts)
    lo = V.min(axis=0); hi = V.max(axis=0)
    ok = False
    for k0 in itertools.product(*[range(-1-int(l), 2-int(h)) for l,h in zip(lo,hi)]):
        if all(tuple(np.array(v)+np.array(k0)) in Aset for v in verts):
            ok = True; break
    if ok:
        tested += 1
        r = in_T([list(g) for g in G])
        found += r
        if tested <= 5 or not r: print('zonotope', G, 'in T:', r, flush=True)
        if tested >= 30: break
print('zonotopes found', tested, 'in T:', found)
