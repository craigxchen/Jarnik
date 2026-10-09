"""Referee (analytic lens) independent checks of Lemma S (lower.md Lemma S = sharp.md Lemma 4.2).

Own Paley-I construction (no code shared with round5/ or sharp_checks/).
  A. H is Hadamard (with the constant column).
  B. Lemma 2.2 (transposed prime-field uncertainty): exhaustive for q = 7, 11 over all ternary
     and {-2..2} zero-sum c with support <= (q-1)/2 (small q), random for q = 43, 47 incl. entries
     divisible by q.
  C. Lemma 2.1 items 2,3,4,5,6 (defect identity, inversion, window, column and half-sum peeling)
     on random integer zero-sum vectors (exact integer arithmetic).
  D. Step-2 small-range bound (2.1): F - M >= f(h) for random ternary c with 4 <= h < 19M/48, and
     the analytic claim min_{4<=h<19M/48, h even} f(h) >= M/16 for every M = q+1, q = 3 mod 4
     prime, 43 <= q <= 10^5.
  E. Lemma S at q = 43: exhaustive over ALL ternary zero-sum c of support 4 and 6 (support 6:
     ~7e7 vectors), plus targeted search near the windows kappa = 1, 2 (columns / half-sums with
     a few entries changed) and hill climbing; reports min F/M over non-exceptional c.
"""
import itertools, math, random, sys
import numpy as np

rng = random.Random(2026)


def chi(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def paley(q):
    """rows: 0..q-1 finite, q = infinity; columns: q labels a = 0..q-1 (nonconstant only)."""
    M = q + 1
    H = np.empty((M, q), dtype=np.int64)
    for x in range(q):
        for a in range(q):
            H[x, a] = -chi(a - x, q) - (1 if a == x else 0)
    H[q, :] = 1
    return H


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


ok = True
out = []

# ---------------- A
for q in (7, 11, 19, 23, 43, 47):
    H = paley(q)
    full = np.hstack([np.ones((q + 1, 1), dtype=np.int64), H])
    good = (full.T @ full == (q + 1) * np.eye(q + 1, dtype=np.int64)).all()
    ok &= bool(good)
out.append(f"A: Hadamard property for q = 7,11,19,23,43,47: {ok}")

# ---------------- B  Lemma 2.2
def check_uncert(q, c):
    H = paley(q)
    T = c @ H
    sp = int(np.count_nonzero(c[:q]))
    if sp > (q - 1) // 2 or not c.any():
        return True, None
    nz = int(np.count_nonzero(T))
    return nz >= (q + 1) // 2 - sp, nz - ((q + 1) // 2 - sp)


worst = {}
for q in (7, 11):
    H = paley(q)
    M = q + 1
    d = (q - 1) // 2
    mn = 10 ** 9
    cnt = 0
    vals = [-2, -1, 1, 2]
    for s in range(1, d + 1):
        for S in itertools.combinations(range(q), s):
            for vv in itertools.product(vals, repeat=s):
                c = np.zeros(M, dtype=np.int64)
                c[list(S)] = vv
                c[q] = -sum(vv)
                T = c @ H
                nz = int(np.count_nonzero(T))
                slack = nz - ((q + 1) // 2 - s)
                mn = min(mn, slack)
                cnt += 1
    worst[q] = (mn, cnt)
    ok &= mn >= 0
for q in (43, 47):
    H = paley(q)
    M = q + 1
    d = (q - 1) // 2
    mn = 10 ** 9
    for it in range(20000):
        s = rng.randint(1, d)
        S = rng.sample(range(q), s)
        c = np.zeros(M, dtype=np.int64)
        for x in S:
            v = rng.choice([-3, -2, -1, 1, 2, 3])
            if it % 3 == 0:
                v *= q
            c[x] = v
        c[q] = -c[:q].sum()
        if it % 5 == 0 and s >= 2:      # force sum over finite part = 0 -> c_inf = 0
            c[S[0]] -= c[:q].sum()
            c[q] = 0
            if c[S[0]] == 0:
                continue
        T = c @ H
        sp = int(np.count_nonzero(c[:q]))
        nz = int(np.count_nonzero(T))
        mn = min(mn, nz - ((q + 1) // 2 - sp))
    worst[q] = (mn, 20000)
    ok &= mn >= 0
out.append(f"B: Lemma 2.2 min slack (#nonzero - ((q+1)/2 - s')): {worst}  (>= 0 required)")

# ---------------- C  Lemma 2.1 items
def F_of(H, c):
    return int(np.abs(c @ H).sum())


cnt = 0
bad = []
for q in (43, 47, 59):
    H = paley(q)
    M = q + 1
    cols = [H[:, a].copy() for a in range(q)]
    for it in range(4000):
        kind = it % 4
        if kind == 0:
            h = rng.randrange(2, M + 1)
            S = rng.sample(range(M), h)
            c = np.zeros(M, dtype=np.int64)
            for x in S:
                c[x] = rng.choice([-3, -2, -1, 1, 2, 3])
            c[S[0]] -= c.sum()
        elif kind == 1:   # perturbed column
            a = rng.randrange(q)
            c = rng.choice([1, -1]) * cols[a].copy()
            for _ in range(rng.randint(1, 3)):
                x, y = rng.sample(range(M), 2)
                c[x] += 1
                c[y] -= 1
        elif kind == 2:   # perturbed half-sum
            a, b2 = rng.sample(range(q), 2)
            c = (rng.choice([1, -1]) * cols[a] + rng.choice([1, -1]) * cols[b2]) // 2
            for _ in range(rng.randint(1, 3)):
                x, y = rng.sample(range(M), 2)
                c[x] += 1
                c[y] -= 1
        else:             # ternary random
            h = 2 * rng.randrange(1, M // 2 + 1)
            S = rng.sample(range(M), h)
            c = np.zeros(M, dtype=np.int64)
            for i, x in enumerate(S):
                c[x] = 1 if i < h // 2 else -1
        if not c.any():
            continue
        assert c.sum() == 0
        T = c @ H
        F = int(np.abs(T).sum())
        n = int(np.abs(c).sum())
        m2 = int((c * c).sum())
        # item 1 parity, Parseval
        ok1 = (T % 2 == 0).all() and int((T * T).sum()) == M * m2 and n % 2 == 0
        # item 2 defect identity (integers): n F = M m2 + sum |T|(n-|T|)
        ok2 = n * F == M * m2 + int((np.abs(T) * (n - np.abs(T))).sum())
        # item 3 inversion
        ok3 = F >= M * int(np.abs(c).max()) and F >= n
        # item 4 window: |F - kappa n| <= 2(F - M)
        kappa = int((2 * np.abs(T) > n).sum())
        ok4 = abs(F - kappa * n) <= 2 * (F - M)
        # item 5 column peeling for every a with c != +-H_a
        ok5 = True
        for a in range(q):
            if (c == cols[a]).all() or (c == -cols[a]).all():
                continue
            if F < min(2 * abs(int(T[a])), 2 * M):
                ok5 = False
        # item 6 half-sum peeling, for the two largest labels
        idx = np.argsort(-np.abs(T))[:2]
        a, b2 = int(idx[0]), int(idx[1])
        s1 = 1 if T[a] >= 0 else -1
        s2 = 1 if T[b2] >= 0 else -1
        g = (s1 * cols[a] + s2 * cols[b2]) // 2
        ok6 = True
        if not (c == g).all():
            phi = lambda t: t - abs(t - M / 2)
            ok6 = F >= M + phi(abs(int(T[a]))) + phi(abs(int(T[b2]))) - 1e-9
        allok = ok1 and ok2 and ok3 and ok4 and ok5 and ok6
        if not allok:
            bad.append((q, it, ok1, ok2, ok3, ok4, ok5, ok6))
        cnt += 1
ok &= not bad
out.append(f"C: Lemma 2.1 items 1-6 on {cnt} vectors (q = 43,47,59): failures = {len(bad)}")

# ---------------- D  small-range bound
def f_small(h, M):
    return (2 * (h - 2) / h) * (M / 2 - h - M / h)


worstD = (10 ** 9, None)
qs = [q for q in range(43, 100001) if q % 4 == 3 and is_prime(q)] if False else None
# faster prime list
N = 100001
sieve = bytearray([1]) * (N + 1)
sieve[0] = sieve[1] = 0
for i in range(2, int(N ** 0.5) + 1):
    if sieve[i]:
        sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
qs = [q for q in range(43, N) if sieve[q] and q % 4 == 3]
for q in qs:
    M = q + 1
    hmax = 19 * M / 48
    best = min(f_small(h, M) for h in range(4, int(math.ceil(hmax)), 2))
    ratio = best / (M / 16)
    if ratio < worstD[0]:
        worstD = (ratio, q)
ok &= worstD[0] >= 1
out.append(f"D1: min over q in [43,1e5] of min_h f(h)/(M/16) = {worstD[0]:.4f} at q = {worstD[1]}  (>= 1 required)")
# also check that (2.1) is a valid inequality on random ternary vectors of small support
viol = 0
tot = 0
for q in (43, 47, 59, 67):
    H = paley(q)
    M = q + 1
    for it in range(3000):
        h = 2 * rng.randrange(2, int(19 * M / 96) + 1)
        if h >= 19 * M / 48:
            continue
        S = rng.sample(range(M), h)
        c = np.zeros(M, dtype=np.int64)
        for i, x in enumerate(S):
            c[x] = 1 if i < h // 2 else -1
        F = F_of(H, c)
        tot += 1
        if F - M < f_small(h, M) - 1e-9:
            viol += 1
ok &= viol == 0
out.append(f"D2: (2.1) F - M >= f(h) on {tot} random small-support ternary vectors: violations = {viol}")

# ---------------- E  exhaustive support 4 and 6 at q = 43
q = 43
H = paley(q)
M = q + 1
minF = {}
for h in (4, 6):
    best = 10 ** 9
    subsets = np.array(list(itertools.combinations(range(M), h)), dtype=np.int32)
    # sign patterns with zero sum, up to global sign: fix first entry +1
    pats = [p for p in itertools.product((1, -1), repeat=h) if sum(p) == 0 and p[0] == 1]
    for p in pats:
        s = np.array(p, dtype=np.int64)
        for start in range(0, len(subsets), 200000):
            Sb = subsets[start:start + 200000]
            T = np.zeros((len(Sb), q), dtype=np.int64)
            for k in range(h):
                T += s[k] * H[Sb[:, k]]
            Fv = np.abs(T).sum(axis=1)
            best = min(best, int(Fv.min()))
    minF[h] = best
ok &= minF[4] >= 17 * M / 16 and minF[6] >= 17 * M / 16
out.append(f"E1: q = 43 exhaustive ternary support 4: min F = {minF[4]}, support 6: min F = {minF[6]}  (17M/16 = {17*M/16:.2f})")

# near-window search: columns and half-sums with small perturbations, plus hill-climb
cols = [H[:, a].copy() for a in range(q)]
def exceptional(c):
    T = c @ H
    aT = np.abs(T)
    if (np.abs(c).sum() == 2):
        return True
    if (aT == M).sum() == 1 and (aT == 0).sum() == q - 1:
        return True
    if (aT == M // 2).sum() == 2 and (aT == 0).sum() == q - 2:
        return True
    return False

best_near = 10 ** 9
for it in range(30000):
    if it % 2 == 0:
        a = rng.randrange(q)
        c = rng.choice([1, -1]) * cols[a].copy()
    else:
        a, b2 = rng.sample(range(q), 2)
        c = (rng.choice([1, -1]) * cols[a] + rng.choice([1, -1]) * cols[b2]) // 2
    for _ in range(rng.randint(1, 4)):
        x, y = rng.sample(range(M), 2)
        c[x] += 1
        c[y] -= 1
    if not c.any() or exceptional(c):
        continue
    best_near = min(best_near, F_of(H, c))
ok &= best_near >= 17 * M / 16
# hill climbing on ternary vectors for low F (non-exceptional)
best_hc = 10 ** 9
for start in range(60):
    h = 2 * rng.randrange(2, M // 2 + 1)
    S = rng.sample(range(M), h)
    c = np.zeros(M, dtype=np.int64)
    for i, x in enumerate(S):
        c[x] = 1 if i < h // 2 else -1
    cur = F_of(H, c) if not exceptional(c) else 10 ** 9
    for step in range(3000):
        c2 = c.copy()
        x, y = rng.sample(range(M), 2)
        c2[x] += 1
        c2[y] -= 1
        if np.abs(c2).max() > 1 or not c2.any() or exceptional(c2):
            continue
        f2 = F_of(H, c2)
        if f2 <= cur:
            c, cur = c2, f2
    best_hc = min(best_hc, cur)
ok &= best_hc >= 17 * M / 16
out.append(f"E2: q = 43 perturbed columns/half-sums min F = {best_near}; ternary hill-climb min F = {best_hc} (non-exceptional; 17M/16 = {17*M/16:.2f})")
out.append("ALL REFEREE LEMMA S CHECKS PASSED" if ok else "FAILURE")
print("\n".join(out))
