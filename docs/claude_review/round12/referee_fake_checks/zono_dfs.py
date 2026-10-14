"""Exhaustive search for the largest lattice zonotope k0 + sum_i [0,L_i] v_i inside a T=1 box (n = 8 models).
DFS over generators in strictly increasing lex order (first nonzero entry > 0), maintaining the feasible
translate set K0 = {k0 : k0 + Z in A}.  Reports zeta(A) (exact, exhaustive) or the best found within a time cap."""
import sys, time, itertools, numpy as np
from rlib import *
name, L = sys.argv[1], int(sys.argv[2]); target = int(sys.argv[3])
n = 12 if name.endswith('11') else 8
pats = {'sdc7': sdc, 'pal7': pal, 'sdc11': sdc, 'pal11': pal}[name](int(name[3:])) if name != 'syl8' else sylvester(3)
A, _ = box(n, pats, L, 0)
d = n - 1
pts = [tuple(x) for x in A.tolist()]
idx = {p: i for i, p in enumerate(pts)}
N = len(pts)
Aarr = np.array(pts)
def shifted(mask, v):
    """mask of k in A with k+v in A and mask[k+v] true"""
    out = np.zeros(N, dtype=bool)
    for i in np.nonzero(mask)[0]:
        q = tuple(Aarr[i] + v)
        j = idx.get(q)
        if j is not None and mask[j]: out[i] = True
    return out
best = [0, None]; t0 = time.time(); nodes = [0]
def canon_pos(v):
    for x in v:
        if x != 0: return x > 0
    return False
def dfs(mask, gens, total, last):
    nodes[0] += 1
    if total > best[0]:
        best[0] = total; best[1] = list(gens)
        print('  total %d with generators %s (%.0fs, nodes %d)' % (total, gens, time.time()-t0, nodes[0]), flush=True)
        if total >= target: raise StopIteration
    sel = Aarr[mask]
    if len(sel) < 2: return
    D = set()
    for a in sel:
        diffs = sel - a
        for v in map(tuple, diffs.tolist()):
            if v > last and canon_pos(v): D.add(v)
    for v in sorted(D):
        va = np.array(v); m1 = mask
        for Lv in range(1, 5):
            m1 = shifted(m1, va * 1) if Lv == 1 else (m1 & shifted_chain(mask, va, Lv))
            if not m1.any(): break
            gens.append((v, Lv)); dfs(m1, gens, total + Lv, v); gens.pop()
def shifted_chain(mask, va, Lv):
    out = mask.copy()
    for j in range(1, Lv + 1):
        out &= shift_by(mask, va * j)
    return out
def shift_by(mask, v):
    out = np.zeros(N, dtype=bool)
    for i in range(N):
        j = idx.get(tuple(Aarr[i] + v))
        if j is not None and mask[j]: out[i] = True
    return out
try:
    dfs(np.ones(N, dtype=bool), [], 0, tuple([-9]*d))
    print('%s L=%d: EXHAUSTIVE, zeta(A) = %d (nodes %d, %.0fs)' % (name, L, best[0], nodes[0], time.time()-t0))
except StopIteration:
    print('%s L=%d: found zonotope of total length %d >= %d: %s (%.0fs)' % (name, L, best[0], target, best[1], time.time()-t0))
