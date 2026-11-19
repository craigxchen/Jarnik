"""Exact counterexample to identifying shortest with degree-four directions."""
from collections import Counter
from itertools import combinations_with_replacement

X = (1,5,8,10,12,16,22,23)
for k in range(1,8):
    moment_signatures = Counter((sum(multiset), sum(x*x for x in multiset))
                                for multiset in combinations_with_replacement(X,k))
    collisions = sum(n*(n-1)//2 for n in moment_signatures.values())
    print(k, collisions)
    assert collisions == (13 if k == 7 else 0)
