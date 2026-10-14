"""Explore the pal11 box (T=1, L=2): size, patterns, symmetry group."""
import itertools, time
import numpy as np
from fplib import pal_patterns, box, canon
M, pats = pal_patterns(11)
print('M', M, 'patterns', len(pats))
for P in pats: print('  ', P)
A = box(M, pats, 2)
print('|A| =', len(A))
Aset = set(A)
d = M - 1
# row permutations fixing base row 0: rows 1+i = i in F_11
def rowperm_from(f):  # f: F_11 -> F_11
    return [0] + [1 + f(i) for i in range(11)]
def act_on_k(perm, k):  # k' indexed by rows 1..11 -> coord x-1
    out = [0]*d
    for x in range(1, M):
        out[perm[x]-1] = k[x-1]
    return tuple(out)
cands = {}
for a in range(1, 11):
    for b in range(11):
        cands[(a, b)] = rowperm_from(lambda i, a=a, b=b: (a*i + b) % 11)
patset = set(pats)
for key, perm in cands.items():
    imgpats = set(canon(M, [perm[x] for x in P]) for P in pats)
    okpat = imgpats == patset
    okA = all(act_on_k(perm, k) in Aset for k in A)
    if okA or okpat:
        print('affine map x->%d x + %d : patterns invariant %s, A invariant %s' % (key[0], key[1], okpat, okA))
# which patterns are active (some k in A with |k(P)| = L)?
for P in pats:
    rel = [x-1 for x in P]
    vals = [sum(k[x] for x in rel) for k in A]
    print(P, 'range of k(P) on A:', min(vals), max(vals))
k0 = [-sum(k) for k in A]
print('base coordinate k_0 range', min(k0), max(k0))
