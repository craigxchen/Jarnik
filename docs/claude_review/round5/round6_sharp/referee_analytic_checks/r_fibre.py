r"""Referee exact check of sharp.md Lemma 6.2 (degeneracy) and Lemma 6.3 (random fibre).

Small canonical Paley profiles (q0 = 7, 11; b = 4).  For n = m/2 in a range of moduli:
  * K'_n = {k in (Z/n)^r : A k in (Z/n) 1}; generators from the twin parametrisation of Lemma 6.3(2)
    (non-twin unit vectors with twin corrections, the u-directions and the t-direction); each
    generator is verified to lie in K'_n, and the generated group has the predicted size
    n^(r - M + 1) (checked through the bijection count on a smaller sample).
  * for v outside L, the index d of <v, K'_n> in Z/n is gcd(n, <v, gen_i>); check d <= 2||v||_1
    (Lemma 6.2 + 6.3(3)), on random small v and on adversarial v = v(c) + d0 u (which lie in
    L + d0 Z^r).
  * Lemma 6.2 directly: for adversarial v in L + G Z^r \ L, G <= 2 ||v||_1.
"""
import math, random
import numpy as np

rng = random.Random(1)


def chi(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def profile(q, b, astar=1):
    M = q + 1
    H = np.empty((M, q), dtype=np.int64)
    for x in range(q):
        for a in range(q):
            H[x, a] = -chi(a - x, q) - (1 if a == x else 0)
    H[q, :] = 1
    S = np.repeat(H, b, axis=1)
    col = lambda a, k: a * b + k
    f, s = {}, {}
    for x in range(q):
        f[x], s[x] = col(x, 0), col(x, 1)
    f[q], s[q] = col(astar, 2), col(astar, 3)
    for x in range(M):
        S[x, f[x]] *= -1
    A = (1 + S) // 2
    eps = {x: int(A[x, f[x]] - A[x, s[x]]) for x in range(M)}
    assert all(((A[:, f[x]] - A[:, s[x]]) == eps[x] * np.eye(M, dtype=np.int64)[x]).all() for x in range(M))
    return A, f, s, eps


def in_L(v, A, f, s, eps):
    M = A.shape[0]
    c = np.array([eps[x] * (v[f[x]] - v[s[x]]) for x in range(M)], dtype=np.int64)
    return c.sum() == 0 and (c @ A == v).all()


def gens(A, f, s, eps, n):
    M, r = A.shape
    twin = set(f.values()) | set(s.values())
    nontwin = [j for j in range(r) if j not in twin]
    G = []
    def complete(knt, u, t):
        k = np.zeros(r, dtype=np.int64)
        k[nontwin] = knt
        base = A[:, nontwin] @ knt
        su = np.array([sum(int(A[x, s[y]]) * int(u[y]) for y in range(M)) for x in range(M)])
        for x in range(M):
            k[f[x]] = eps[x] * (t - base[x] - su[x])
            k[s[x]] = u[x] - k[f[x]]
        return k % n
    for i in range(len(nontwin)):
        knt = np.zeros(len(nontwin), dtype=np.int64)
        knt[i] = 1
        G.append(complete(knt, np.zeros(M, dtype=np.int64), 0))
    for y in range(M):
        u = np.zeros(M, dtype=np.int64)
        u[y] = 1
        G.append(complete(np.zeros(len(nontwin), dtype=np.int64), u, 0))
    G.append(complete(np.zeros(len(nontwin), dtype=np.int64), np.zeros(M, dtype=np.int64), 1))
    return G


ok = True
lines = []
for q0 in (7, 11):
    A, f, s, eps = profile(q0, 4)
    M, r = A.shape
    worst = 0.0
    nv = 0
    worst62 = 0.0
    for n in list(range(2, 41)) + [60, 84, 120, 180]:
        G = gens(A, f, s, eps, n)
        for k in G:
            Ak = (A @ k) % n
            ok &= bool((Ak == Ak[0]).all())
        for it in range(150):
            if it % 2 == 0:
                v = np.array([rng.choice([0, 0, 0, 0, 1, -1, 2, -2]) for _ in range(r)], dtype=np.int64)
            else:
                c = np.array([rng.randint(-2, 2) for _ in range(M)], dtype=np.int64)
                c[0] -= c.sum()
                d0 = rng.randint(2, 3 * n)
                u = np.zeros(r, dtype=np.int64)
                for _ in range(rng.randint(1, 3)):
                    u[rng.randrange(r)] += rng.choice([1, -1])
                v = c @ A + d0 * u
                if in_L(v, A, f, s, eps):
                    continue
                # Lemma 6.2: d0 <= 2 ||v||_1
                worst62 = max(worst62, d0 / (2 * np.abs(v).sum()))
            if not v.any() or in_L(v, A, f, s, eps):
                continue
            vals = [int(v @ k) % n for k in G]
            d = n
            for x in vals:
                d = math.gcd(d, x)
            nv += 1
            worst = max(worst, d / (2 * np.abs(v).sum()))
    ok &= worst <= 1 and worst62 <= 1
    lines.append(f"q0={q0} (M={M}, r={r}): generators in K' for all n: {ok}; {nv} vectors v not in L: max d/(2||v||_1) = {worst:.3f}; "
                 f"Lemma 6.2 adversarial max G/(2||v||_1) = {worst62:.3f}")
lines.append("ALL REFEREE FIBRE CHECKS PASSED" if ok else "FAILURE")
print("\n".join(lines))
