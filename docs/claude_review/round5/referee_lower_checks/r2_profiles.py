"""Referee check r2: profiles of lower.md Sec. 3 (exact integers).

Canonical assignment (Prop 3.2) and random one-flip capacity-two assignments (Prop 3.3).
Built from the text: label a has b physical copies of H_a; row x flipped in one copy of label alpha(x);
distinct rows use distinct physical columns; each flipped column gets a distinct unflipped sibling.
Checks:
  * identity (3.1) 2B = bF + sum_x(|T_alpha(x) - 2 c_x H(x,alpha(x))| - |T_alpha(x)|) against direct c^T S;
  * twins S_f(x) - S_s(x) = -2 H(x,alpha(x)) e_x, rank S = M, every column nonconstant;
  * Prop 3.2 (canonical, b = 5, 8): pairs 2B-r >= b-4, ternary support h: 2B-r >= b+2h-8,
    amplitude A>=2: 2B >= (b-2) M A; exhaustive pairs and support 4 (q = 43, b = 8), random otherwise;
    reports whether the bounds are attained;
  * Prop 3.3 (random capacity-two, b = 8, 35): table rows, and the intermediate 2B - r >= b(F-M+1) - 2n.
  * twisted extension (referee): for sum c != 0, ||S^T c||_1 >= max(2n, (b-2)(M-n)).
"""
import random, itertools
import numpy as np
from rl_common import paley_H, T_of, Fval, in_E

random.seed(5)
np.random.seed(5)
FAIL = []


def build(H, b, alpha):
    """alpha: list of labels (0..q-1) per row; capacity <= 2; returns S, flipcol, sibcol"""
    M = H.shape[0]; q = M - 1
    assert max(alpha.count(a) for a in set(alpha)) <= 2 and b >= 4
    cols = []; flipcol = {}; sibcol = {}
    for a in range(q):
        rows = [x for x in range(M) if alpha[x] == a]
        base = len(cols)
        for t in range(b):
            col = H[:, a + 1].copy()
            if t < len(rows):
                col[rows[t]] *= -1
                flipcol[rows[t]] = base + t
            cols.append(col)
        for t, x in enumerate(rows):
            sibcol[x] = base + len(rows) + t   # distinct unflipped siblings
    S = np.array(cols, dtype=np.int64).T
    return S, flipcol, sibcol


def twoB(S, c):
    return int(np.abs(c @ S).sum())


def identity31(H, S, b, alpha, c):
    T = T_of(H, c)
    F = int(np.abs(T).sum())
    flip = sum(abs(int(T[alpha[x]]) - 2 * int(c[x]) * int(H[x, alpha[x] + 1])) - abs(int(T[alpha[x]]))
               for x in range(len(c)))
    return b * F + flip


def structure(H, S, flipcol, sibcol, alpha):
    M = H.shape[0]
    for x in range(M):
        d = S[:, flipcol[x]] - S[:, sibcol[x]]
        e = np.zeros(M, dtype=np.int64); e[x] = -2 * H[x, alpha[x] + 1]
        if not (d == e).all():
            FAIL.append(("twin", x))
    if np.linalg.matrix_rank(S.astype(float)) != M:
        FAIL.append(("rank",))
    if any(len(set(S[:, j])) == 1 for j in range(S.shape[1])):
        FAIL.append(("constant column",))
    if len(set(map(tuple, S))) != M:
        FAIL.append(("rows not distinct",))


def ternary_random(M, h):
    supp = random.sample(range(M), h)
    c = np.zeros(M, dtype=np.int64); c[supp[:h // 2]] = 1; c[supp[h // 2:]] = -1
    return c


report = []
# ---------------- canonical
for q, b in ((11, 5), (11, 8), (19, 8), (43, 8), (43, 5), (47, 8)):
    H = paley_H(q); M = q + 1; r = b * q
    astar = 0
    alpha = list(range(q)) + [astar]
    S, fc, sc = build(H, b, alpha)
    structure(H, S, fc, sc, alpha)
    # identity on random vectors
    for _ in range(400):
        c = np.random.randint(-3, 4, size=M).astype(np.int64); c[-1] -= c.sum()
        if twoB(S, c) != identity31(H, S, b, alpha, c):
            FAIL.append(("3.1", q, b))
    mins = {"pair": 10 ** 9, "tern": 10 ** 9, "amp": 10 ** 9}
    for x, y in itertools.combinations(range(M), 2):
        c = np.zeros(M, dtype=np.int64); c[x] = 1; c[y] = -1
        mins["pair"] = min(mins["pair"], twoB(S, c) - r - (b - 4))
    # exhaustive support 4 via row sums (direct c^T S)
    if q <= 47:
        for x, y in itertools.combinations(range(M), 2):
            P = S[x] + S[y]
            rest = [z for z in range(M) if z not in (x, y)]
            for z, w in itertools.combinations(rest, 2):
                v = int(np.abs(P - S[z] - S[w]).sum())
                mins["tern"] = min(mins["tern"], v - r - (b + 2 * 4 - 8))
    for _ in range(3000):
        h = random.randrange(4, M + 1, 2)
        c = ternary_random(M, h)
        mins["tern"] = min(mins["tern"], twoB(S, c) - r - (b + 2 * h - 8))
    for a in range(1, M):
        c = H[:, a].copy()
        mins["tern"] = min(mins["tern"], twoB(S, c) - r - (b + 2 * M - 8))
    for a, a2 in itertools.combinations(range(1, M), 2):
        for s in (1, -1):
            c = (H[:, a] + s * H[:, a2]) // 2
            h = int((c != 0).sum())
            mins["tern"] = min(mins["tern"], twoB(S, c) - r - (b + 2 * h - 8))
    for _ in range(1500):
        A = random.randint(2, 4)
        c = np.random.randint(-A, A + 1, size=M).astype(np.int64)
        c[0] = A
        c[-1] -= c.sum()
        if c[-1] > A or c[-1] < -A:
            continue
        A = int(np.abs(c).max())
        mins["amp"] = min(mins["amp"], twoB(S, c) - (b - 2) * M * A)
    for k, v in mins.items():
        if v < 0:
            FAIL.append(("Prop3.2", q, b, k, v))
    report.append(("canonical", q, b, mins))

# ---------------- random capacity-two assignments
for q, b in ((43, 8), (47, 8), (43, 35), (59, 8)):
    H = paley_H(q); M = q + 1; r = b * q
    for rep in range(3):
        slots = [a for a in range(q) for _ in range(2)]
        random.shuffle(slots)
        alpha = slots[:M]
        S, fc, sc = build(H, b, alpha)
        structure(H, S, fc, sc, alpha)
        mins = {"pair": 10 ** 9, "col": 10 ** 9, "half": 10 ** 9, "other": 10 ** 9, "amp": 10 ** 9,
                "inter": 10 ** 9}
        for x, y in itertools.combinations(range(M), 2):
            c = np.zeros(M, dtype=np.int64); c[x] = 1; c[y] = -1
            mins["pair"] = min(mins["pair"], twoB(S, c) - r - (b - 4))
        for a in range(1, M):
            for s in (1, -1):
                c = s * H[:, a]
                mins["col"] = min(mins["col"], twoB(S, c) - r - (b + 2 * M - 8))
        for a, a2 in itertools.combinations(range(1, M), 2):
            for s in (1, -1):
                for t in (1, -1):
                    c = (s * H[:, a] + t * H[:, a2]) // 2
                    mins["half"] = min(mins["half"], twoB(S, c) - r - (b + M - 16))
        for _ in range(3000):
            h = random.randrange(4, M + 1, 2)
            c = ternary_random(M, h)
            if in_E(H, c):
                continue
            v = twoB(S, c) - r
            F = Fval(H, c)
            mins["inter"] = min(mins["inter"], v - (b * (F - M + 1) - 2 * h))
            mins["other"] = min(mins["other"], v - (b * (M / 16 + 1) - 2 * h))
        for _ in range(800):
            A = 2
            c = np.random.randint(-1, 2, size=M).astype(np.int64)
            c[0] = 2; c[-1] -= c.sum()
            if abs(c[-1]) > 2:
                continue
            A = int(np.abs(c).max())
            mins["amp"] = min(mins["amp"], twoB(S, c) - r - (b * (A * M - M + 1) - 2 * A * M))
        for k, v in mins.items():
            if v < 0:
                FAIL.append(("Prop3.3", q, b, k, v))
        report.append(("cap2", q, b, rep, mins))

# ---------------- twisted extension (sum c != 0)
for q, b in ((11, 8), (19, 8), (43, 8)):
    H = paley_H(q); M = q + 1
    alpha = list(range(q)) + [0]
    S, fc, sc = build(H, b, alpha)
    worst = 10 ** 9
    for _ in range(4000):
        n_target = random.randint(1, 2 * M)
        c = np.zeros(M, dtype=np.int64)
        for _ in range(n_target):
            c[random.randrange(M)] += random.choice([1, -1])
        if c.sum() == 0 or c.sum() % 2 or not c.any():
            continue
        n = int(np.abs(c).sum())
        v = twoB(S, c)
        worst = min(worst, v - max(2 * n, (b - 2) * (M - n)))
    if worst < 0:
        FAIL.append(("twisted", q, worst))
    report.append(("twisted ||S^T c||_1 - max(2n,(b-2)(M-n))", q, b, worst))

for line in report:
    print(line)
print("ALL OK" if not FAIL else "FAILURES: %s" % FAIL[:10])
