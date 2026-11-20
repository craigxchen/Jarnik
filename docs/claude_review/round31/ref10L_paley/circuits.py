"""For Paley S (order m), deleted row set M (|M|=k), T = [1; S_M] (rank k+1).
Count r x r column minors of T that vanish (r=k+1): if none vanish, every kernel vector of T has
support >= r+1, so the divisor method of X.3 cannot use smaller column sets."""
from itertools import combinations
import sys
sys.path.insert(0, '.')
from paley_checks import paley_S, det_int, rank_Q
for m in (11, 19):
    S = paley_S(m)
    for k in range(0, (m+1)//2):
        if m == 19 and k > 3: break
        sets = singular = 0; minsupp_reduced = 0
        for M in combinations(range(m), k):
            T = [[1]*m] + [S[i] for i in M]
            r = k + 1
            assert rank_Q(T) == r
            sets += 1
            z = sum(1 for c in combinations(range(m), r)
                    if det_int([[row[j] for j in c] for row in T]) == 0)
            singular += z
            if z: minsupp_reduced += 1
        print(f"m={m} k={k}: {sets} row sets; {singular} vanishing {k+1}x{k+1} column minors in total; "
              f"{minsupp_reduced} row sets admit a kernel vector of support <= {k+1}")
        sys.stdout.flush()
