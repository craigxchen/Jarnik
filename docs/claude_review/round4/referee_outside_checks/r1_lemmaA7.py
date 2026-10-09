"""Referee check of Lemma A.7 (outside.md) on an INDEPENDENT build of the fake-circle profile.

Differences from the author's check_fake_barrier.py:
  * flip rows, designated flips and siblings are assigned at random (subject to the stated rules),
    not cyclically;
  * the weighted inequality (iii) is tested with ACTUAL primes = 1 mod 4 in [P, 2P] whose O_a log-sums
    are balanced into a window of length 1/M (by swap hill-climbing), logs at 50 digits (decimal);
  * besides exhaustive small characters and random ones, an adversarial local search tries to
    minimise the slack of (iii) and of its equal-weight form.
Checked statements (n = ||c||_1, c zero-sum integral, nonzero):
  (i)   F(c) = sum_{a>=1} |T_a| >= M
  (ii)  twin structure: designated flipped column minus sibling = -2 H(x,a) e_x  (=> saturation), rank S = M
  (iii) V_c - W/2 >= tau (n + 0.86 M),  V_c = (1/2) sum_j w_j |c^T S_j|,  W = sum_j w_j,  tau = log P
  (iii_eq) equal weights: sum_j |c^T S_j| - r >= 2 (n + M)
"""
import random
import sys
from decimal import Decimal, getcontext
from itertools import combinations, product
import numpy as np

getcontext().prec = 50
random.seed(2026)


def sylvester(m):
    H = np.array([[1]], dtype=np.int64)
    while H.shape[0] < m:
        H = np.block([[H, H], [H, -H]])
    return H


def legendre(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1


def paley(q):
    n = q + 1
    S = np.zeros((n, n), dtype=np.int64)
    S[0, 1:] = 1
    S[1:, 0] = -1
    for i in range(q):
        for j in range(q):
            S[1 + i, 1 + j] = legendre(j - i, q)
    return np.eye(n, dtype=np.int64) + S


def normalize(H):
    H = H * H[:, [0]]
    M = H.shape[0]
    assert (H.T @ H == M * np.eye(M, dtype=np.int64)).all() and (H[:, 0] == 1).all()
    return H


def build(H, b):
    """columns: label a=1..M-1, copies 0..b-1; F_a = copies {0,1}, S_a = {2,3}, O_a = rest."""
    M = H.shape[0]
    labels = list(range(1, M))
    # flipped columns: (a,0),(a,1); assign rows: every row 1 or 2 flips, random
    flipped = [(a, t) for a in labels for t in (0, 1)]
    random.shuffle(flipped)
    rows = list(range(M)) + random.sample(range(M), len(flipped) - M)  # each row once, M-2 rows twice
    random.shuffle(rows)
    flip_row = dict(zip(flipped, rows))
    # designated flip per row (injective): choose one flipped column per row
    designated = {}
    for col, x in flip_row.items():
        if x not in designated or random.random() < 0.5:
            designated[x] = col
    assert len(designated) == M and len(set(designated.values())) == M
    # siblings: bijection F_a -> S_a, random
    sib = {}
    for a in labels:
        perm = [2, 3]
        random.shuffle(perm)
        sib[(a, 0)] = (a, perm[0])
        sib[(a, 1)] = (a, perm[1])
    cols = [(a, t) for a in labels for t in range(b)]
    idx = {col: j for j, col in enumerate(cols)}
    S = np.zeros((M, len(cols)), dtype=np.int64)
    for j, (a, t) in enumerate(cols):
        S[:, j] = H[:, a]
        if (a, t) in flip_row:
            S[flip_row[(a, t)], j] *= -1
    fr = [list(flip_row.values()).count(x) for x in range(M)]
    assert min(fr) >= 1 and max(fr) <= 2
    # twin structure
    for x, col in designated.items():
        d = S[:, idx[col]] - S[:, idx[sib[col]]]
        e = np.zeros(M, dtype=np.int64)
        e[x] = -2 * H[x, col[0]]
        assert (d == e).all()
    assert np.linalg.matrix_rank(S.astype(float)) == M
    return S, cols, idx


def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def weights(cols, M, b, P):
    """Assign distinct primes = 1 mod 4 in [P,2P] (random sample) to columns; balance O_a sums."""
    r = len(cols)
    pool = []
    n = P + (1 - P) % 4
    # sample candidate primes spread over [P,2P]
    while len(pool) < 3 * r:
        m = random.randrange(P, 2 * P)
        m += (1 - m) % 4
        if m < 2 * P and is_prime(m) and m not in pool:
            pool.append(m)
    chosen = random.sample(pool, r)
    logs = {p: Decimal(p).ln() for p in pool}
    labels = list(range(1, M))
    k = b - 4
    # O_a sets: start with random partition of the first (M-1)k chosen primes
    Oprimes = chosen[:(M - 1) * k]
    rest = chosen[(M - 1) * k:]
    O = {a: Oprimes[(a - 1) * k:a * k] for a in labels}
    spare = [p for p in pool if p not in chosen]

    def sums():
        return {a: sum(logs[p] for p in O[a]) for a in labels}
    target = Decimal(1) / M
    s = sums()
    it = 0
    while max(s.values()) - min(s.values()) > target:
        it += 1
        amax = max(s, key=s.get)
        amin = min(s, key=s.get)
        gap = s[amax] - s[amin]
        best = None
        # swap between amax and amin, or swap amax element with spare prime
        for i, p in enumerate(O[amax]):
            for j, q in enumerate(O[amin]):
                d = logs[p] - logs[q]
                if d > 0:
                    newgap = abs(gap - 2 * d)
                    if best is None or newgap < best[0]:
                        best = (newgap, 'pair', i, j)
            for j, q in enumerate(spare[:200]):
                d = logs[p] - logs[q]
                if 0 < d < gap:
                    newgap = max(abs(gap - d), Decimal(0))
                    if best is None or newgap < best[0]:
                        best = (newgap, 'spare', i, j)
        if best is None or best[0] >= gap:
            random.shuffle(spare)
            if it > 5000:
                raise RuntimeError("balancing failed")
            continue
        _, kind, i, j = best
        if kind == 'pair':
            O[amax][i], O[amin][j] = O[amin][j], O[amax][i]
        else:
            O[amax][i], spare[j] = spare[j], O[amax][i]
        s = sums()
    # assign
    wcol = [None] * r
    ri = iter(rest)
    for j, (a, t) in enumerate(cols):
        if t >= 4:
            wcol[j] = logs[O[a][t - 4]]
        else:
            wcol[j] = logs[next(ri)]
    used = [p for a in labels for p in O[a]] + rest
    assert len(set(used)) == r
    return wcol, max(s.values()) - min(s.values())


def characters(M, maxn, nrand):
    cs = []
    for k in range(2, maxn + 1):
        for supp in combinations(range(M), k):
            for vals in product(*[range(-maxn, maxn + 1)] * k):
                if 0 in vals or sum(vals) != 0 or sum(abs(v) for v in vals) > maxn:
                    continue
                c = np.zeros(M, dtype=np.int64)
                c[list(supp)] = vals
                cs.append(c)
    for _ in range(nrand):
        c = np.zeros(M, dtype=np.int64)
        for _ in range(random.randint(1, 3 * M)):
            i, j = random.sample(range(M), 2)
            amp = random.choice([1, 1, 1, 2, 3])
            c[i] += amp
            c[j] -= amp
        if c.any():
            cs.append(c)
    return cs


def run(H, b, name, P, maxn, nrand, nsearch):
    M = H.shape[0]
    S, cols, idx = build(H, b)
    r = S.shape[1]
    w, spread = weights(cols, M, b, P)
    tau = Decimal(P).ln()
    W = sum(w)
    wf = np.array([float(x) for x in w])
    cs = characters(M, maxn, nrand)
    for a in range(1, M):
        cs.append(H[:, a].copy())
        for a2 in range(a + 1, M):
            for s in (1, -1):
                cs.append((H[:, a] + s * H[:, a2]) // 2)
    wf = np.array([float(x) for x in w])
    worst_w = None
    worst_e = None
    for c in cs:
        n = int(np.abs(c).sum())
        T = c @ H
        assert T[0] == 0
        F = int(np.abs(T[1:]).sum())
        assert F >= M
        col = np.abs(c @ S)
        Geq = int(col.sum()) - r
        se = Geq / 2 - (n + M)
        assert se >= 0, (name, c, Geq)
        swf = (0.5 * float(col @ wf) - float(W) / 2) / float(tau) - (n + 0.86 * M)
        if swf < 1e-6:  # re-verify in 50-digit arithmetic
            V = sum(Decimal(int(col[j])) * w[j] for j in range(r) if col[j]) / 2
            swf = (V - W / 2) / tau - (n + Decimal("0.86") * M)
        sw = swf
        assert sw >= 0, (name, c, sw)
        worst_e = se if worst_e is None else min(worst_e, se)
        worst_w = sw if worst_w is None else min(worst_w, sw)
    # adversarial local search on the weighted slack (float, then re-verified in Decimal)
    def slack_f(c):
        n = np.abs(c).sum()
        col = np.abs(c @ S)
        V = 0.5 * float(col @ wf)
        return (V - float(W) / 2) / float(tau) - (n + 0.86 * M)
    best_found = None
    for _ in range(nsearch):
        c = np.zeros(M, dtype=np.int64)
        i, j = random.sample(range(M), 2)
        c[i], c[j] = 1, -1
        cur = slack_f(c)
        for _ in range(300):
            i, j = random.sample(range(M), 2)
            c2 = c.copy()
            c2[i] += 1
            c2[j] -= 1
            if not c2.any():
                continue
            v2 = slack_f(c2)
            if v2 <= cur or random.random() < 0.05:
                c, cur = c2, v2
        if best_found is None or cur < best_found[0]:
            best_found = (cur, c.copy())
    c = best_found[1]
    n = int(np.abs(c).sum())
    col = np.abs(c @ S)
    V = sum(Decimal(int(col[j])) * w[j] for j in range(r) if col[j]) / 2
    sw = (V - W / 2) / tau - (n + Decimal("0.86") * M)
    assert sw >= 0
    print("%-9s M=%2d b=%3d r=%4d P=%d tau=%.2f window=%.4f(<=1/M=%.4f): %6d chars OK; "
          "min weighted slack %.3f, min equal-weight slack %.1f; local search min slack %.3f (n=%d)"
          % (name, M, b, r, P, float(tau), float(spread), 1 / M, len(cs), float(worst_w), worst_e,
             float(sw), n))
    sys.stdout.flush()


def main():
    run(normalize(sylvester(8)), 48, "Sylvester", 10 ** 6, 6, 3000, 300)
    run(normalize(paley(11)), 72, "Paley", 10 ** 6, 4, 3000, 300)
    run(normalize(sylvester(16)), 96, "Sylvester", 10 ** 7, 4, 2000, 200)
    run(normalize(paley(19)), 120, "Paley", 10 ** 7, 4, 1000, 100)


if __name__ == "__main__":
    main()
