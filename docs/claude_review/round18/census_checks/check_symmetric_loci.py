"""Prop S check: for every k in 5..8 and every label permutation sigma != id of {0..k-1}, there is a
boundary split T (2 <= |T| <= k-2) such that T and sigma(T) are incompatible (so D_T and D_{sigma T}
are disjoint).  Also report, per conjugacy class (cycle type), the number of boundary divisors
D_T with T incompatible to sigma(T) (these have d_T bounded on the fixed locus Y_sigma)."""
from itertools import combinations, permutations
def splits(k):
    out = set()
    U = frozenset(range(k))
    for s in range(2, k - 1):
        for S in combinations(range(k), s):
            S = frozenset(S)
            out.add(frozenset([S, U - S]))
    return out
def compatible(A, B, U):
    a = next(iter(A)); b = next(iter(B))
    ac = U - a; bc = U - b
    return (not (a & b)) or (not (a & bc)) or (not (ac & b)) or (not (ac & bc))
for k in range(5, 9):
    U = frozenset(range(k))
    sp = list(splits(k))
    bytype = {}
    worst = None
    for sig in permutations(range(k)):
        if all(sig[i] == i for i in range(k)):
            continue
        ctype = []
        seen = set()
        for i in range(k):
            if i in seen: continue
            c = 0; j = i
            while j not in seen:
                seen.add(j); j = sig[j]; c += 1
            ctype.append(c)
        ctype = tuple(sorted(ctype, reverse=True))
        if ctype in bytype:
            continue
        cnt = 0
        for s in sp:
            img = frozenset(frozenset(sig[x] for x in part) for part in s)
            if img != s and not compatible(s, img, U):
                cnt += 1
        bytype[ctype] = cnt
        assert cnt > 0, (k, sig)
    print("k=%d  #splits=%d  min #(T incompatible with sigma T) over sigma!=id = %d" % (k, len(sp), min(bytype.values())))
    print("   by cycle type:", dict(sorted(bytype.items())))
