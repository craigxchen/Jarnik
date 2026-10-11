"""Referee check r6 (exact integers): is Lemma S needed for Theorem 5.1(b)?

Research-notes item 393 (paley_primefield_polynomial_character_gap.md) already proves, for EVERY one-flip
capacity-two assignment with b >= 5 copies, two bounds on 2B - r (B = sum_j |c^T S_j|/2, r = b(M-1)) for a
zero-sum ternary c of support h:
  (i)  [polynomial/uncertainty, h < q]  F >= h(2M - 2h)/(h + 2)  (from #{T_a != 0} >= M/2 - h and
       (t-2)(h-t) >= 0), hence  2B - r >= b(F - M + 1) - 2h >= b(h(2M-2h)/(h+2) - M + 1) - 2h;
  (ii) [high support, h >= 8]  2B - r >= b + 2h - 8M/h     (k <= 2 #{T_a != 0}, capacity two).
Referee consequence: for M >= 70 every non-pair ternary c has 2B - r >= 0.6 M - 22 (b >= 5), with no
exceptional set, versus Lemma S's (b/16 - 2) M + b, which is positive only for b > 32.
This script checks (i), (ii) and min over h of max((i),(ii)) on random capacity-two assignments
(q = 43, 47, 59, 83; b = 5, 8), exhaustively at h = 4 for q = 43, on random ternary c of every even
support, on all columns and half-sums, and reports the observed minimum of (2B - r)/M over non-pair
ternary characters.
"""
import random, itertools
import numpy as np
from rl_common import paley_H, T_of

random.seed(393)
np.random.seed(393)
FAIL = []


def build(H, b, alpha):
    M = H.shape[0]; q = M - 1
    cols = []
    for a in range(q):
        rows = [x for x in range(M) if alpha[x] == a]
        for t in range(b):
            col = H[:, a + 1].copy()
            if t < len(rows):
                col[rows[t]] *= -1
            cols.append(col)
    return np.array(cols, dtype=np.int64).T


def bound_i(M, b, h):
    # valid for h < q: F >= h(2M-2h)/(h+2); 2B - r >= b(F - M + 1) - 2h  (rational lower bound)
    F_lb = h * (2 * M - 2 * h) / (h + 2)
    return b * (F_lb - M + 1) - 2 * h


def bound_ii(M, b, h):
    return b + 2 * h - 8 * M / h if h >= 8 else -10 ** 9


report = []
for q in (43, 47, 59, 83):
    H = paley_H(q); M = q + 1
    for b in (5, 8):
        r = b * q
        slots = [a for a in range(q) for _ in range(2)]
        random.shuffle(slots)
        alpha = slots[:M]
        S = build(H, b, alpha)
        worst_i = worst_ii = 10 ** 9
        worst_ratio = 10 ** 9
        worst_comb = 10 ** 9
        def test(c):
            global worst_i, worst_ii, worst_ratio, worst_comb
            h = int((c != 0).sum())
            v = int(np.abs(c @ S).sum()) - r
            T = T_of(H, c)
            F = int(np.abs(T).sum())
            s = int((T != 0).sum())
            if h < q:
                if s < M // 2 - h:
                    FAIL.append(("uncertainty", q, h))
                if F * (h + 2) < h * (2 * M - 2 * h):
                    FAIL.append(("F_i", q, h))
                bi = bound_i(M, b, h)
                worst_i = min(worst_i, v - bi)
            if h >= 8:
                bii = bound_ii(M, b, h)
                worst_ii = min(worst_ii, v - bii)
            comb = max(bound_i(M, b, h) if h < q else -10 ** 9, bound_ii(M, b, h))
            worst_comb = min(worst_comb, v - comb)
            worst_ratio = min(worst_ratio, v / M)
        if q == 43:
            for x, y in itertools.combinations(range(M), 2):
                for z, w in itertools.combinations([t for t in range(M) if t not in (x, y)], 2):
                    if (z, w) < (x, y):
                        continue
                    c = np.zeros(M, dtype=np.int64); c[[x, y]] = 1; c[[z, w]] = -1
                    test(c)
        for _ in range(6000):
            h = random.randrange(4, M + 1, 2)
            supp = random.sample(range(M), h)
            c = np.zeros(M, dtype=np.int64); c[supp[:h // 2]] = 1; c[supp[h // 2:]] = -1
            test(c)
        for a in range(1, M):
            test(H[:, a].copy())
        for a, a2 in itertools.combinations(range(1, M), 2):
            for t in (1, -1):
                test((H[:, a] + t * H[:, a2]) // 2)
        # analytic minimum over h of the combined bound
        comb_min = min(max(bound_i(M, b, h) if h < q else -10 ** 9, bound_ii(M, b, h))
                       for h in range(4, M + 1, 2))
        lemmaS = (b / 16 - 2) * M + b
        if worst_i < 0 or worst_ii < 0 or worst_comb < 0:
            FAIL.append(("bound", q, b, worst_i, worst_ii, worst_comb))
        report.append("q=%d M=%d b=%d: slack over (i) >= %s, over (ii) >= %s; min_h max((i),(ii)) = %.1f "
                      "(= %.2f M); Lemma-S row (b/16-2)M+b = %.1f; observed min (2B-r)/M over non-pair "
                      "ternary = %.3f" % (q, M, b, worst_i, worst_ii, comb_min, comb_min / M, lemmaS,
                                          worst_ratio))
for line in report:
    print(line)
# asymptotic form of the combined bound for b = 5 at large M (exact rational evaluation)
for M in (10 ** 3, 10 ** 4, 10 ** 5):
    b = 5
    cm = min(max(bound_i(M, b, h), bound_ii(M, b, h)) for h in range(4, M - 1, 2))
    print("b=5, M=%d: min_h max((i),(ii)) / M = %.4f" % (M, cm / M))
print("ALL OK" if not FAIL else "FAILURES %s" % FAIL[:5])
