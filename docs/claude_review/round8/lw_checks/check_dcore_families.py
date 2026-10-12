"""Structured characters of the Paley skew double core, and its Liouville margins.

Part A: s(c) for (i) columns of K (c = K(.,a), zero-sum for a >= 1), (ii) label pairs
        (K(.,a) - K(.,b))/2 and (iii) row pairs e_x - e_y (first 12 indices).  q = 7, 11, 19, 23, 43, 47.
Part B: Liouville margin of exact.md Thm 4.1(iii): Phi(c) = sum_j |<c,s_j>| >= 2(M-1) + ||c||_1 on
        random integer zero-sum c (entries up to 6, all supports) -- an independent spot check.
Part C: support-weight fraction w_S/W = s/(#types) at the minimum, compared with 1/4.
Part D: the K^T family equals the K family with rows 0 and a flipped (exact).
Usage: python3 check_dcore_families.py
"""
import random
import numpy as np
from dcore import paley_skew, double_core_types

random.seed(7)
out = []
def log(s):
    print(s); out.append(s)

for q in (7, 11, 19, 23, 43, 47):
    K = paley_skew(q); M = q + 1
    f1, f2 = double_core_types(K)
    T = np.array(f1 + f2, dtype=np.int64)
    ntypes = T.shape[0]
    def s(c):
        return int(np.count_nonzero(T @ np.array(c, dtype=np.int64)))
    cols = [[K[x][a] for x in range(M)] for a in range(1, M)]
    s_cols = sorted(set(s(c) for c in cols))
    lp = []
    for a in range(1, min(M, 12)):
        for b in range(a + 1, min(M, 12)):
            c = [(K[x][a] - K[x][b]) // 2 for x in range(M)]
            lp.append(s(c))
    rp = []
    for x in range(min(M, 12)):
        for y in range(x + 1, min(M, 12)):
            c = [0] * M; c[x] = 1; c[y] = -1
            rp.append(s(c))
    log(f"q={q} M={M} types={ntypes}: s(columns of K) = {s_cols} (M/2+1 = {M//2+1}); "
        f"label pairs: min {min(lp)} max {max(lp)}; row pairs: min {min(rp)} max {max(rp)}; "
        f"min-fraction s/#types = {min(s_cols+lp+rp)/ntypes:.4f} vs 1/4")
    # Part B
    worst = None
    for trial in range(3000 if q < 40 else 600):
        k = random.randint(2, M)
        supp = random.sample(range(M), k)
        c = [0] * M
        for x in supp[:-1]:
            c[x] = random.choice([-1, 1]) * random.randint(1, 6)
        c[supp[-1]] = -sum(c)
        if all(v == 0 for v in c):
            continue
        phi = int(np.abs(T @ np.array(c, dtype=np.int64)).sum())
        slack = phi - (2 * (M - 1) + sum(abs(v) for v in c))
        if worst is None or slack < worst:
            worst = slack
        assert slack >= 0, (q, c)
    log(f"   Liouville margin spot check: min over random c of Phi - (2(M-1)+||c||_1) = {worst} (>= 0)")

# Part D: the second family is the first family modified on rows {0, a} (section 5.1 of lw.md)
for q in (7, 11, 19, 23, 43, 47):
    K = paley_skew(q); M = q + 1
    f1, f2 = double_core_types(K)
    ok = all(f2[a-1][x] == (-f1[a-1][x] if x in (0, a) else f1[a-1][x]) for a in range(1, M) for x in range(M))
    assert ok
    log(f"Part D q={q}: cut S'_a = cut S_a with rows 0 and a flipped, for every a: {ok}")

open('check_dcore_families_out.txt', 'w').write('\n'.join(out) + '\n')
