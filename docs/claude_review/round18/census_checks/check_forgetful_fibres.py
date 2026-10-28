"""Cor. 7.3 count: for k = 8 and every Q subset {0..7} with 4 <= |Q| <= 7, the number of boundary
splits T of {0..7} whose restriction to Q is a boundary split of Q (both sides >= 2)."""
from itertools import combinations
k = 8
U = frozenset(range(k))
splits = set()
for s in range(2, k - 1):
    for S in combinations(range(k), s):
        S = frozenset(S); splits.add(frozenset([S, U - S]))
for q in range(4, 8):
    counts = set()
    for Q in combinations(range(k), q):
        Q = frozenset(Q)
        c = 0
        for sp in splits:
            a, b = tuple(sp)
            if len(a & Q) >= 2 and len(b & Q) >= 2:
                c += 1
        counts.add(c)
    print("|Q| = %d: #splits restricting to a boundary split of Q: %s" % (q, sorted(counts)))
