"""Vertical-line templates of degree 8: z_sigma(u) = prod_{k=1..8} (u + i sigma_k b_k), b_k distinct positive.
Bounded C (families.md Lemma P, D = 8: Im p_r agree for r = 1, 2, 3) iff sum sigma_k b_k and
sum sigma_k b_k^3 are constant (even r give real power sums).  Search 8-sets b in {1..B}:
fibres of (subset sum, subset cube sum) with >= 6 members, and their affine rank."""
import itertools, sys
import numpy as np
B = int(sys.argv[1]) if len(sys.argv) > 1 else 20
D = 8
best = {}
for b in itertools.combinations(range(1, B + 1), D):
    fib = {}
    for mask in range(1 << D):
        s1 = s3 = 0
        for k in range(D):
            if mask >> k & 1:
                s1 += b[k]; s3 += b[k] ** 3
        fib.setdefault((s1, s3), []).append(mask)
    for key, masks in fib.items():
        if len(masks) >= 4:
            vecs = np.array([[1 if (m >> k) & 1 else -1 for k in range(D)] for m in masks], dtype=float)
            rk = int(np.linalg.matrix_rank(vecs[1:] - vecs[0]))
            kk = (len(masks), rk)
            if kk not in best:
                best[kk] = (b, key)
for kk in sorted(best):
    print("fibre size %d, affine rank %d: example b=%s, (sum, cube sum) = %s" % (kk[0], kk[1], best[kk][0], best[kk][1]))
