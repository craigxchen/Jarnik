r"""Exhaustive test in the complete Fano block Y_3(2) (t = 3, j = 2), degree e = 7 = 3.5 j:
which products of e heavy2 restricted generators (coset differences with |H'| = 2, entries +-2 on two axes) lie in
J(Y_3(2))^perp_7 ?  All multisets with axis counts <= 2j = 4 (necessary, grid) and line counts N_line <= nmax are tested.
mod-p negatives are rigorous (T_e mod p contains the reduction of the rational T_e).
usage: python3 yblock_exh.py e nmax
"""
import sys, itertools, time
import numpy as np
from math import factorial
from flint import nmod_mat
sys.argv_saved = list(sys.argv)
import yblock as yb
PR = yb.PR
e, nmax = int(sys.argv[1]), int(sys.argv[2])
t, j = 3, 2
M, Y = yb.build_Y(t, j)
b = 7
Zs = list(range(8))
heavy = []
seen = set()
for h in range(1, 8):
    Hs = [0, h]
    cosets = sorted({tuple(sorted({x ^ y for y in Hs})) for x in Zs})
    for C1, C2 in itertools.combinations(cosets, 2):
        c = np.zeros(8, dtype=np.int64)
        for x in C1: c[x] += 1
        for x in C2: c[x] -= 1
        a = tuple(M @ c)
        if a not in seen and tuple(-x for x in a) not in seen:
            seen.add(a); heavy.append(np.array(a))
lines = [(a, c, a ^ c) for a in range(1, 8) for c in range(a + 1, 8) if (a ^ c) > c]
# tangent space T_e mod p
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
    for ex, cf in poly.items():
        vec[idx[ex]] = (cf % PR) * minv[idx[ex]] % PR
    for i in range(rT):
        c = int(vec[piv[i]])
        if c: vec = (vec - (c * Ta[i]) % PR) % PR
    return not vec.any()
nz = [np.nonzero(g)[0] for g in heavy]
tl = [[li for li, L in enumerate(lines) if any(heavy[i][q - 1] for q in L)] for i in range(len(heavy))]
cnt = np.zeros(b, dtype=np.int64); lcnt = np.zeros(len(lines), dtype=np.int64)
stats = {}; tested = [0]; members = []
t0 = time.time()
def rec(start, ch):
    if len(ch) == e:
        tested[0] += 1
        mem = member([heavy[i] for i in ch])
        key = (int(cnt.max()), int(lcnt.max()))
        a, tt = stats.get(key, (0, 0)); stats[key] = (a + int(mem), tt + 1)
        if mem and len(members) < 10: members.append(tuple(ch))
        return
    for i in range(start, len(heavy)):
        cnt[nz[i]] += 1; lcnt[tl[i]] += 1
        if (cnt <= 2 * j).all() and (lcnt <= nmax).all():
            ch.append(i); rec(i, ch); ch.pop()
        cnt[nz[i]] -= 1; lcnt[tl[i]] -= 1
rec(0, [])
print(f'Y_3(2), e={e}, dim T_e (mod p)={rT}; heavy2 multisets with axis<=4, lines<={nmax}: tested {tested[0]}  ({time.time()-t0:.0f}s)')
print('(max axis count, max line count) -> (in J^perp, total):')
for k, v in sorted(stats.items()): print('   ', k, v)
for mm in members[:5]: print('  member example:', [list(heavy[i]) for i in mm])
