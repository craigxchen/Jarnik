"""Independent exact check of the canonical-assignment bounds (research notes item 392,
paley_canonical_all_character_half_height.md), re-proved in lower.md Sec. 3.2:

Canonical profile: Paley H (convention of paley.py), b copies per nonconstant label, finite row x
flipped at one copy of label x (the label whose transform coordinate is T_x), row infinity flipped at
a second copy of a fixed label a*.  B(c) = sum_j |c^T S_j|/2, r = b(M-1).
  pairs (h = 2):                 B - r/2 >= b/2 - 2
  ternary, support h >= 4:       B - r/2 >= b/2 + h - 4
  amplitude A = max|c_x| >= 2:   B >= M A (b/2 - 1)
NOTE: in paley.py the diagonal entry is H(x, label x) = -chi(0) - 1 = -1, as in the note."""
import random
from itertools import combinations, product
import numpy as np
from paley import paley

random.seed(11)
FAIL = []


def build(q, b, astar=0):
    H = paley(q)
    M = q + 1
    cols = []
    for a in range(q):          # label a in F_q is column index a+1
        flips = [a]             # finite row a flipped in a copy of label a
        if a == astar:
            flips.append(q)     # infinity row flipped in a second copy of label a*
        for t in range(b):
            col = H[:, a + 1].copy()
            if t < len(flips):
                col[flips[t]] = -col[flips[t]]
            cols.append(col)
    S = np.array(cols, dtype=np.int64).T
    assert S.shape == (M, b * q)
    return H, S


def B2(S, c):  # 2B
    return int(np.abs(c @ S).sum())


def run(q, b):
    H, S = build(q, b)
    M = q + 1
    r = b * q
    worst = {"pair": 10 ** 9, "tern": 10 ** 9, "amp": 10 ** 9}
    # pairs, exhaustive
    for x, y in combinations(range(M), 2):
        c = np.zeros(M, dtype=np.int64); c[x] = 1; c[y] = -1
        d = B2(S, c) - r          # = 2(B - r/2)
        worst["pair"] = min(worst["pair"], d - (b - 4))
    # ternary: exhaustive h=4 (and h=6 for q<=11), random larger
    tested = 0
    hs = [4, 6] if q <= 11 else [4]
    for h in hs:
        for supp in combinations(range(M), h):
            for vals in product([1, -1], repeat=h):
                if sum(vals) != 0 or vals[0] < 0:
                    continue
                c = np.zeros(M, dtype=np.int64); c[list(supp)] = vals
                d = B2(S, c) - r
                worst["tern"] = min(worst["tern"], d - (b + 2 * h - 8))
                tested += 1
    for _ in range(20000):
        h = random.randrange(4, M + 1, 2)
        supp = random.sample(range(M), h)
        c = np.zeros(M, dtype=np.int64); c[supp[: h // 2]] = 1; c[supp[h // 2:]] = -1
        d = B2(S, c) - r
        worst["tern"] = min(worst["tern"], d - (b + 2 * h - 8))
        tested += 1
    # columns and half-sums explicitly (ternary, h = M and M/2)
    for a in range(1, M):
        c = H[:, a].copy()
        d = B2(S, c) - r
        worst["tern"] = min(worst["tern"], d - (b + 2 * M - 8))
    for a, a2 in combinations(range(1, M), 2):
        for s in (1, -1):
            c = (H[:, a] + s * H[:, a2]) // 2
            h = int((c != 0).sum())
            d = B2(S, c) - r
            worst["tern"] = min(worst["tern"], d - (b + 2 * h - 8))
    # amplitude >= 2
    for _ in range(5000):
        s = random.randint(2, M)
        supp = random.sample(range(M), s)
        vals = [random.randint(-4, 4) for _ in range(s - 1)]
        vals.append(-sum(vals))
        c = np.zeros(M, dtype=np.int64); c[supp] = vals
        A = int(np.abs(c).max())
        if A < 2:
            continue
        worst["amp"] = min(worst["amp"], B2(S, c) - M * A * (b - 2))   # 2B - 2MA(b/2-1)
    for k, v in worst.items():
        if v < 0:
            FAIL.append((q, b, k, v))
    print("q=%3d b=%d: min margins (>=0 required): pairs %d, ternary %d (%d tested), amplitude %d"
          % (q, b, worst["pair"], worst["tern"], tested, worst["amp"]))


if __name__ == "__main__":
    for q, b in ((11, 5), (19, 5), (23, 8), (43, 5), (47, 8), (59, 5)):
        run(q, b)
    print("ALL CHECKS PASSED" if not FAIL else "FAILURES: %s" % FAIL)
