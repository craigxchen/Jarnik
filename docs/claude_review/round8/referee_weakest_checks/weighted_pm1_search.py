"""Referee check: can {0,+-1} vertical templates refute the WEIGHTED repair SLW_w(k, g < 1/4) for k >= 6?

A {0,+-1} template with D columns and s odd-power conditions (sum sigma_l b_l^(2j-1) = 0, j <= s) has
weighted exponent gamma_w -> (2s+1)/(2D).  For k = 6 and gamma_w >= 1/4 one needs s >= 2 and D <= 4s+2,
with 6 independent {0,+-1} solutions.  This script computes, for s = 2 and D = 8, 9, 10 and every
b in binom({1..B}, D), the rank of the set of {0,+-1} solutions of sum sigma b = sum sigma b^3 = 0.
Usage: python3 weighted_pm1_search.py D B
"""
import sys
import itertools
import numpy as np


def main(D, B, batch=128):
    V = np.array([v for v in itertools.product((-1, 0, 1), repeat=D) if any(v)], dtype=np.int64)
    # keep one of each +-pair
    first = np.array([v[np.flatnonzero(v)[0]] for v in V])
    V = V[first > 0]
    best = (0, None)
    hist = {}
    combos = itertools.combinations(range(1, B + 1), D)
    total = 0
    while True:
        chunk = list(itertools.islice(combos, batch))
        if not chunk:
            break
        bs = np.array(chunk, dtype=np.int64)
        s1 = V @ bs.T
        s3 = V @ (bs ** 3).T
        Z = (s1 == 0) & (s3 == 0)
        for j in range(bs.shape[0]):
            idx = np.flatnonzero(Z[:, j])
            rk = 0 if len(idx) == 0 else int(np.linalg.matrix_rank(V[idx].astype(float)))
            hist[rk] = hist.get(rk, 0) + 1
            if rk > best[0]:
                best = (rk, chunk[j], V[idx].tolist())
        total += len(chunk)
    print(f"D={D}, s=2, b subsets of {{1..{B}}} ({total} sets): histogram of rank of {{0,+-1}} solutions: "
          f"{dict(sorted(hist.items()))}")
    print(f"  max rank {best[0]} at b={best[1]}; gamma_w of a rank-r template = 5/(2D) = {5/(2*D):.4f}")
    if best[0] and len(best[2]) <= 12:
        for v in best[2]:
            print("    ", v)


if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))
