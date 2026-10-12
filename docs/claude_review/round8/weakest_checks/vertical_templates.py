"""Vertical-line polynomial templates z_sigma(u) = prod_k (u + i sigma_k b_k), b_k distinct positive.
Bounded C (families.md Lemma P, D = 6) iff sum_k sigma_k b_k is constant on the family.
Search: 6-sets b in {1..B} with >= 6 equal-sum subsets; report the affine rank of the sign vectors."""
import itertools, sys
import numpy as np
B = int(sys.argv[1]) if len(sys.argv) > 1 else 30
D = int(sys.argv[2]) if len(sys.argv) > 2 else 6
best = {}
for b in itertools.combinations(range(1, B + 1), D):
    sums = {}
    for mask in range(1 << D):
        s = sum(b[k] for k in range(D) if mask >> k & 1)
        sums.setdefault(s, []).append(mask)
    for s, masks in sums.items():
        if len(masks) >= 6:
            vecs = np.array([[1 if (m >> k) & 1 else -1 for k in range(D)] for m in masks], dtype=float)
            rk = int(np.linalg.matrix_rank(vecs[1:] - vecs[0]))
            key = (len(masks), rk)
            if key not in best:
                best[key] = (b, s)
for key in sorted(best):
    print("fibre size %d, affine rank %d: example b=%s, subset sum %d" % (key[0], key[1], best[key][0], best[key][1]))
