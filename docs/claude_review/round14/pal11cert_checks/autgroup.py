"""Full automorphism group of the pal11 box A (L=2) among signed permutations of the 11 coordinates of k'.
Method: backtracking; an automorphism must preserve A; we check A-invariance at the end (exact set test) and
prune with the invariant 'pattern incidence' (sets P with |k(P)|<=2 forced).  Exhaustive."""
import itertools, sys
import numpy as np
from fplib import pal_patterns, box
M, pats = pal_patterns(11)
A = box(M, pats, 2); Aset = set(A); d = 11
Aarr = np.array(A)
# invariant: for subsets T of coordinates of size <= 3 and sign vectors, the number of k in A with given restriction
# cheap invariant: projection counts onto pairs (x,y): N[(x,sx),(y,sy)] = #{k in A: k_x = sx, k_y = sy}
cnt = {}
for x in range(d):
    for y in range(d):
        if x == y: continue
        for sx in (-1, 1):
            for sy in (-1, 1):
                cnt[(x, sx, y, sy)] = int(((Aarr[:, x] == sx) & (Aarr[:, y] == sy)).sum())
autos = []
def rec(img):  # img: list of (coord, sign) images of coords 0..len-1
    i = len(img)
    if i == d:
        perm = img
        ok = True
        for k in A:
            kk = [0]*d
            for x in range(d):
                c, s = perm[x]; kk[c] = s * k[x]
            if tuple(kk) not in Aset: ok = False; break
        if ok: autos.append(tuple(perm))
        return
    used = {c for c, s in img}
    for c in range(d):
        if c in used: continue
        for s in (1, -1):
            good = True
            for x in range(i):
                cx, sx_ = img[x]
                for a in (-1, 1):
                    for b in (-1, 1):
                        if cnt[(x, a, i, b)] != cnt[(cx, sx_*a, c, s*b)]: good = False; break
                    if not good: break
                if not good: break
            if good:
                img.append((c, s)); rec(img); img.pop()
rec([])
print('number of signed-permutation automorphisms of A:', len(autos))
for g in autos: print(g)
