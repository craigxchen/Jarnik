# invariant factors of the anchor-difference lattice of actual best clusters
from aux_checks import *
from collections import Counter
cnt = Counter(); rk = Counter()
for circ in CIRCLES:
    for m in (3,4,5,6,7,8):
        for C, zs, avecs, units in clusters(circ, m, 4):
            B = [[a - b for a, b in zip(avecs[i], avecs[0])] for i in range(1, m)]
            U, S, V = smith(B)
            d = [S[i][i] for i in range(min(len(S), len(S[0]))) if S[i][i] != 0]
            cnt[(m, max(d))] += 1; rk[(m, len(d))] += 1
print("largest invariant factor counts (m, d_max):", sorted(cnt.items()))
print("rank counts (m, rank):", sorted(rk.items()))
