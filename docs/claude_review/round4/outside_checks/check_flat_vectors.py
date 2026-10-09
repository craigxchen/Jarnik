"""Exact census of 'flat' characters (F(c) = ||H^T c||_1 = M) for Paley-I and Sylvester matrices:
all zero-sum c with ||c||_1 <= 6 and entries in [-3,3], plus columns, two-column half-sums and
random sparse combinations of columns.  Reports flat vectors that are NOT pairs, columns or
two-column half-sums (supports the remark in outside.md Sec. 2.2; evidence only for larger
supports)."""
import random
from itertools import combinations, product
import numpy as np
from check_fake_barrier import sylvester, paley1, normalize

random.seed(9)


def census(H, name, maxn=6):
    M = H.shape[0]
    seen = 0
    other = []
    def F(c):
        T = c @ H
        assert T[0] == 0
        return int(np.abs(T[1:]).sum())
    for k in range(2, maxn + 1):
        for supp in combinations(range(M), k):
            for vals in product([-3, -2, -1, 1, 2, 3], repeat=k):
                if sum(vals) != 0 or sum(abs(v) for v in vals) > maxn:
                    continue
                c = np.zeros(M, dtype=np.int64)
                c[list(supp)] = vals
                seen += 1
                if F(c) == M and k > 2:
                    other.append(c)
        if seen > 400000:
            break
    # combos of few columns
    for _ in range(4000):
        k = random.randint(1, 4)
        cols = random.sample(range(1, M), k)
        signs = [random.choice([1, -1]) for _ in cols]
        v = sum(s * H[:, a] for s, a in zip(signs, cols))
        for d in (1, 2, 4):
            if (v % d == 0).all():
                c = v // d
                if c.any() and F(c) == M:
                    nz = int((c != 0).sum())
                    if nz not in (M, M // 2):
                        other.append(c)
    print("%-10s M=%2d: %d small characters scanned; flat non-pair vectors of support <= %d: %d; "
          "other flat column-combinations: %d"
          % (name, M, seen, maxn, sum(1 for c in other if (c != 0).sum() <= maxn),
             sum(1 for c in other if (c != 0).sum() > maxn)))


for q in (11, 19, 23):
    census(normalize(paley1(q)), "Paley-I")
for m in (8, 16):
    census(normalize(sylvester(m)), "Sylvester")
