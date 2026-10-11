"""Referee check r1: Lemma 2.1, Lemma 2.2 and Lemma S of lower.md (exact integers).

(A) skew structure K = D H = I + J, J skew (used in Prop 3.2).
(B) Lemma 2.2 (transposed uncertainty): #{a: T_a != 0} >= (q+1)/2 - s' for s' <= (q-1)/2,
    exhaustive at q = 7 for all c with entries in [-2,2] and finite support <= 3, random elsewhere,
    including adversarial supports (quadratic-residue sets, arithmetic progressions) and q-divisible c.
(C) Lemma 2.1 (2)-(6) on random c (defect identity, integrality, peelings).
(D) Lemma S adversarial search at q = 43, 47, 59: exhaustive support 4 (all zero-sum ternary),
    random ternary of every even support, structured starts (perturbed columns / half-sums /
    triple and quadruple combinations, rounded projections) followed by greedy descent on F over
    non-exceptional ternary vectors with moves preserving zero sum.
    Reports min (F-M)/M over non-exceptional c and compares with 1/16.
"""
import sys, random, itertools
import numpy as np
from rl_common import paley_H, T_of, Fval, in_E, chi

random.seed(20261009)
np.random.seed(7)
FAIL = []

# ---------------- (A) skew structure
for q in (7, 11, 19, 43, 47):
    H = paley_H(q)
    M = q + 1
    D = np.diag([-1] * q + [1])
    K = D @ H  # rows: x finite get -1, infinity +1 ; columns: const, labels
    # reorder columns so that the constant column is identified with infinity (index q)
    perm = list(range(1, M)) + [0]
    Kp = K[:, perm]
    J = Kp - np.eye(M, dtype=np.int64)
    ok = (J.T == -J).all() and (np.abs(J) + np.eye(M, dtype=np.int64) == 1).all()
    if not ok:
        FAIL.append(("skew", q))
print("(A) K = DH = I + J with J skew, |J_ij| = 1 off diagonal: q = 7,11,19,43,47:",
      "ok" if not FAIL else FAIL)

# ---------------- (B) Lemma 2.2
def nonzero_count(H, c):
    return int((T_of(H, c) != 0).sum())

worst = {}
cnt = 0
q = 7
H = paley_H(q); M = q + 1
for s in range(1, 4):
    for supp in itertools.combinations(range(q), s):
        for vals in itertools.product([v for v in range(-2, 3) if v != 0], repeat=s):
            for cinf in range(-6, 7):
                c = np.zeros(M, dtype=np.int64)
                c[list(supp)] = vals
                c[q] = cinf
                if c.sum() != 0:
                    continue
                cnt += 1
                gap = nonzero_count(H, c) - ((q + 1) // 2 - s)
                worst[q] = min(worst.get(q, 99), gap)
                if gap < 0:
                    FAIL.append(("L2.2", q, c.tolist()))
for q in (11, 19, 43, 47, 59):
    H = paley_H(q); M = q + 1
    d = (q - 1) // 2
    QR = [x for x in range(1, q) if chi(x, q) == 1]
    for trial in range(4000):
        mode = trial % 4
        s = random.randint(1, d)
        if mode == 0:
            supp = random.sample(range(q), s)
        elif mode == 1:
            supp = random.sample(QR, min(s, len(QR)))
        elif mode == 2:
            st = random.randrange(q); step = random.randrange(1, q)
            supp = sorted({(st + k * step) % q for k in range(s)})
        else:
            supp = random.sample(range(q), s)
        s = len(supp)
        c = np.zeros(M, dtype=np.int64)
        for x in supp:
            v = 0
            while v == 0:
                v = random.randint(-3, 3)
            c[x] = v * (q if mode == 3 else 1)
        c[q] = -int(c[:q].sum())
        cnt += 1
        gap = nonzero_count(H, c) - ((q + 1) // 2 - s)
        worst[q] = min(worst.get(q, 99), gap)
        if gap < 0:
            FAIL.append(("L2.2", q, c.tolist()))
print("(B) Lemma 2.2: %d vectors; min slack per q:" % cnt, worst,
      "FAIL" if any(f[0] == "L2.2" for f in FAIL) else "ok")

# ---------------- (C) Lemma 2.1
def check21(H, c):
    M = H.shape[0]
    T = T_of(H, c)
    F = int(np.abs(T).sum()); n = int(np.abs(c).sum()); m2 = int((c * c).sum())
    assert (T % 2 == 0).all() and n % 2 == 0 and (np.abs(T) <= n).all()
    assert int((T * T).sum()) == M * m2
    # defect identity, times n
    assert n * F == M * m2 + int((np.abs(T) * (n - np.abs(T))).sum())
    assert F >= M and F >= n and F >= M * int(np.abs(c).max())
    kappa = int((2 * np.abs(T) > n).sum())
    assert abs(F - kappa * n) <= 2 * (F - M)
    is_col = any((c == H[:, a]).all() or (c == -H[:, a]).all() for a in range(1, M))
    for a in range(M - 1):
        ta = abs(int(T[a]))
        if not ((c == H[:, a + 1]).all() or (c == -H[:, a + 1]).all()):
            assert F >= min(2 * ta, 2 * M)
    # half-sum peeling at the two largest |T|
    idx = np.argsort(-np.abs(T))[:2]
    a, b = idx[0], idx[1]
    s = 1 if T[a] >= 0 else -1; t = 1 if T[b] >= 0 else -1
    g = (s * H[:, a + 1] + t * H[:, b + 1]) // 2
    if not (c == g).all():
        phi = lambda x: x - abs(x - M / 2)
        assert F >= M + phi(abs(int(T[a]))) + phi(abs(int(T[b]))) - 1e-9
    return True

n21 = 0
for q in (11, 19, 43, 47):
    H = paley_H(q); M = q + 1
    for trial in range(3000):
        kind = trial % 3
        if kind == 0:
            h = random.randrange(2, M + 1, 2)
            supp = random.sample(range(M), h)
            c = np.zeros(M, dtype=np.int64); c[supp[:h // 2]] = 1; c[supp[h // 2:]] = -1
        elif kind == 1:
            c = np.random.randint(-2, 3, size=M).astype(np.int64)
            c[-1] -= c.sum()
        else:
            a, b = random.sample(range(1, M), 2)
            c = (H[:, a] + random.choice([1, -1]) * H[:, b]) // 2
            x, y = random.sample(range(M), 2)
            c = c.copy(); c[x] += 1; c[y] -= 1
        if not c.any():
            continue
        check21(H, c); n21 += 1
print("(C) Lemma 2.1 (1)-(6): %d vectors ok" % n21)

# ---------------- (D) Lemma S search
def ternary_moves(c):
    """zero-sum preserving local moves on ternary vectors"""
    M = len(c)
    x, y = random.sample(range(M), 2)
    d = c.copy()
    d[x] += 1; d[y] -= 1
    if np.abs(d).max() > 1:
        return None
    return d

def greedy(H, c, steps):
    M = H.shape[0]
    best = c.copy(); bF = Fval(H, best)
    for _ in range(steps):
        d = ternary_moves(best)
        if d is None or not d.any() or in_E(H, d):
            continue
        Fd = Fval(H, d)
        if Fd <= bF:
            best, bF = d, Fd
    return best, bF

summary = {}
for q in (43, 47, 59):
    H = paley_H(q); M = q + 1
    best = (10 ** 9, None)
    # exhaustive support 4: c = e_x + e_y - e_z - e_w
    Hrows = H[:, 1:]
    pairsum = {}
    nsupp4 = 0
    mn4 = 10 ** 9
    for x, y in itertools.combinations(range(M), 2):
        P = Hrows[x] + Hrows[y]
        for z, w in itertools.combinations(range(M), 2):
            if len({x, y, z, w}) < 4 or (z, w) < (x, y):
                continue
            T = P - Hrows[z] - Hrows[w]
            F = int(np.abs(T).sum()); nsupp4 += 1
            mn4 = min(mn4, F)
    best = min(best, (mn4, "supp4"))
    # random ternary, every even support
    for trial in range(6000):
        h = random.randrange(6, M + 1, 2)
        supp = random.sample(range(M), h)
        c = np.zeros(M, dtype=np.int64); c[supp[:h // 2]] = 1; c[supp[h // 2:]] = -1
        if in_E(H, c):
            continue
        F = Fval(H, c)
        if F < best[0]:
            best = (F, "random h=%d" % h)
    # structured starts + greedy descent
    starts = []
    for _ in range(60):
        a, b, cc, dd = random.sample(range(1, M), 4)
        s = [random.choice([1, -1]) for _ in range(4)]
        v1 = H[:, a].copy()                              # column, perturbed
        k = random.randint(1, 4)
        for _ in range(k):
            x, y = random.sample(range(M), 2)
            v1[x] = 0 if v1[x] != 0 else v1[x]; v1[y] = 0 if v1[y] != 0 else v1[y]
        v1[-1] -= v1.sum() if abs(v1[-1] - v1.sum()) <= 1 else 0
        starts.append(v1)
        v2 = (s[0] * H[:, a] + s[1] * H[:, b]) // 2       # half-sum, perturbed
        x, y = random.sample(range(M), 2)
        v2 = v2.copy(); v2[x] = 0; v2[y] = 0
        starts.append(v2)
        tri = s[0] * H[:, a] + s[1] * H[:, b] + s[2] * H[:, cc]   # sign of a 3-combination
        starts.append(np.sign(tri) * (np.abs(tri) == 3))
        quad = s[0] * H[:, a] + s[1] * H[:, b] + s[2] * H[:, cc] + s[3] * H[:, dd]
        starts.append(np.sign(quad) * (np.abs(quad) == 4))
        starts.append(np.sign(quad) * (np.abs(quad) >= 2))
    nst = 0
    for c in starts:
        c = np.clip(c, -1, 1).astype(np.int64)
        # repair zero sum
        while c.sum() != 0:
            sgn = 1 if c.sum() > 0 else -1
            idx = [i for i in range(M) if c[i] == sgn]
            c[random.choice(idx)] = 0
        if not c.any():
            continue
        c, F = greedy(H, c, 1500)
        nst += 1
        if not in_E(H, c) and c.any():
            if F < best[0]:
                best = (F, "greedy |c|_1=%d" % int(np.abs(c).sum()))
    ratio = (best[0] - M) / M
    summary[q] = (best[0], round(ratio, 4), best[1], nsupp4, nst)
    if ratio < 1 / 16:
        FAIL.append(("LemmaS", q, best))
print("(D) Lemma S search: q -> (min F over non-exceptional, (F-M)/M, where, #supp4, #greedy runs)")
for q, v in summary.items():
    print("    q=%d M=%d:" % (q, q + 1), v)
print("ALL OK" if not FAIL else ("FAILURES: %s" % FAIL[:5]))
