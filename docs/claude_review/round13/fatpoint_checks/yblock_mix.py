r"""Mixed light/heavy products in J(Y_3(2))^perp_e: test the 'weighted (zonotope) cost' rule
   sum_i w_i [axis a touched by generator i] <= 2j  (w = 1 light, 2 heavy2)
on random multisets with plain axis counts <= 2j and line counts <= 3j (necessary).  Negatives rigorous (mod p).
usage: python3 yblock_mix.py e nsamp
"""
import sys, itertools, numpy as np
from math import factorial
from flint import nmod_mat
import yblock as yb
PR = yb.PR
e, nsamp = int(sys.argv[1]), int(sys.argv[2])
t, j, b = 3, 2, 7
M, Y = yb.build_Y(t, j)
Zs = list(range(8)); gens = []; wts = []; seen = set()
for z1, z2 in itertools.combinations(Zs, 2):
    c = np.zeros(8, dtype=np.int64); c[z1], c[z2] = 1, -1
    gens.append(M @ c); wts.append(1)
for h in range(1, 8):
    cos = sorted({tuple(sorted({x ^ y for y in [0, h]})) for x in Zs})
    for C1, C2 in itertools.combinations(cos, 2):
        c = np.zeros(8, dtype=np.int64)
        for x in C1: c[x] += 1
        for x in C2: c[x] -= 1
        a = tuple(M @ c)
        if a not in seen and tuple(-x for x in a) not in seen:
            seen.add(a); gens.append(np.array(a)); wts.append(2)
lines = [(a, c, a ^ c) for a in range(1, 8) for c in range(a + 1, 8) if (a ^ c) > c]
lowm = []
for f in range(e): lowm += yb.monos(b, f)
X, nul = yb.to_nmod(yb.evalm(Y, lowm)).nullspace()
me = yb.monos(b, e)
Em = yb.to_nmod(yb.evalm(Y, me).T)
Vm = nmod_mat(nul, X.nrows(), [int(X[i, k]) for k in range(nul) for i in range(X.nrows())], PR)
T = (Vm * Em).rref()[0]; rT = T.rank()
Ta = np.array([[int(T[i, c]) for c in range(T.ncols())] for i in range(rT)], dtype=np.int64)
piv = [int(np.nonzero(Ta[i])[0][0]) for i in range(rT)]
idx = {a: i for i, a in enumerate(me)}
minv = [pow(factorial(e) // int(np.prod([factorial(int(x)) for x in a])) % PR, -1, PR) for a in me]
def member(G):
    poly = {tuple([0] * b): 1}
    for g in G:
        new = {}
        for ex, cf in poly.items():
            for x in range(b):
                if g[x]:
                    e2 = list(ex); e2[x] += 1; e2 = tuple(e2)
                    new[e2] = new.get(e2, 0) + cf * int(g[x])
        poly = {k: v for k, v in new.items() if v}
    vec = np.zeros(len(me), dtype=np.int64)
    for ex, cf in poly.items(): vec[idx[ex]] = (cf % PR) * minv[idx[ex]] % PR
    for i in range(rT):
        c = int(vec[piv[i]])
        if c: vec = (vec - (c * Ta[i]) % PR) % PR
    return not vec.any()
rng = np.random.default_rng(7)
stats = {}
nz = [np.nonzero(g)[0] for g in gens]
tl = [[li for li, L in enumerate(lines) if any(gens[i][q - 1] for q in L)] for i in range(len(gens))]
for s in range(nsamp):
    nh = int(rng.integers(0, e + 1))  # number of heavy
    for _ in range(50):
        ch = list(rng.integers(28, len(gens), size=nh)) + list(rng.integers(0, 28, size=e - nh))
        cnt = np.zeros(b, int); wc = np.zeros(b, int); lc = np.zeros(7, int)
        for i in ch: cnt[nz[i]] += 1; wc[nz[i]] += wts[i]; lc[tl[i]] += 1
        if (cnt <= 2 * j).all() and (lc <= 3 * j).all(): break
    else:
        continue
    mem = member([gens[i] for i in ch])
    key = (nh, int(wc.max()) <= 2 * j)
    a, tt = stats.get(key, (0, 0)); stats[key] = (a + int(mem), tt + 1)
print(f'Y_3(2) e={e} dimT={rT}: (number of heavy, weighted cost <= 2j) -> (in J^perp, sampled)')
for k, v in sorted(stats.items()): print('   ', k, v)
