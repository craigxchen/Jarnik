"""Exact checks for Lemma S (l^1-stability of the Paley transform) and its ingredients.

Part 1  Uncertainty U^T (transposed form of round4/flipped.md Lemma 6.2, with the infinity row):
        for c in Z^M, sum c = 0, s' = |supp c cap F_q| <= (q-1)/2:
        #{a in F_q : T_a != 0} >= (q+1)/2 - s'.          exhaustive (q=7,11) + random.
Part 2  Identities/inequalities used in the proof, on random and structured c:
        (D) defect identity  F = M m2/n + (1/n) sum |T_a|(n-|T_a|);
        (I) integrality      |F - k n| <= 2(F - M),  k = #{a : |T_a| > n/2};
        (L) F >= n;
        (P1) F >= min(2|T_a|, 2M) for every label a with c != +-H_a;
        (P2) F >= M + phi(T_a) + phi(sigma T_b) for c != (H_a + sigma H_b)/2 (T_a>0, sigma=sgn T_b),
             phi(t) = |t| - |t - M/2|;
        (S) small range      F - M >= (2(n-2)/n)(M/2 - n - M/n) for 4 <= n <= (q-1)/2 (n = ||c||_1).
Part 3  Search for near-flat non-exceptional characters: exhaustive small supports, structured
        families, hill climbing.  Lemma S claims (F - M)/M >= 1/16 for q >= 43.
All F, T are exact integers."""
import random
import sys
from itertools import combinations, product
import numpy as np
from paley import paley, transform, F, exceptional, is_prime

random.seed(20261009)
np.random.seed(20261009)
FAIL = []


def check(cond, msg):
    if not cond:
        FAIL.append(msg)
        print("FAIL:", msg)


def nonzero_count(H, c):
    return int((transform(H, c) != 0).sum())


# ---------------------------------------------------------------- Part 1
def part1():
    print("== Part 1: transposed prime-field uncertainty ==")
    tested = 0
    for q in (7, 11):
        H = paley(q)
        M = q + 1
        d = (q - 1) // 2
        rng = range(-2, 3) if q == 7 else range(-1, 2)
        maxs = 8 if q == 7 else 6
        minslack = 10 ** 9
        for s in range(1, maxs + 1):
            for supp in combinations(range(M), s):
                for vals in product([v for v in rng if v != 0], repeat=s):
                    if sum(vals) != 0:
                        continue
                    c = np.zeros(M, dtype=np.int64)
                    c[list(supp)] = vals
                    sp = sum(1 for x in supp if x < q)
                    if sp > d:
                        continue
                    nz = nonzero_count(H, c)
                    tested += 1
                    check(nz >= (q + 1) // 2 - sp, "U^T q=%d c=%s" % (q, c))
                    minslack = min(minslack, nz + sp - (q + 1) // 2)
        print("q=%d exhaustive: min (#nonzero + s' - (q+1)/2) = %d" % (q, minslack))
    for q in (19, 23, 31, 43, 47, 59, 67, 71, 79, 83):
        H = paley(q)
        M = q + 1
        d = (q - 1) // 2
        minslack = 10 ** 9
        for _ in range(3000):
            s = random.randint(2, d)
            supp = random.sample(range(M), s)
            vals = [random.choice([1, -1, 2, -2, 3, q, -q, 7 * q * q]) for _ in range(s - 1)]
            vals.append(-sum(vals))
            if vals[-1] == 0:
                continue
            c = np.zeros(M, dtype=np.int64)
            c[supp] = vals
            sp = sum(1 for x in supp if x < q)
            if sp > d:
                continue
            nz = nonzero_count(H, c)
            tested += 1
            check(nz >= (q + 1) // 2 - sp, "U^T q=%d" % q)
            minslack = min(minslack, nz + sp - (q + 1) // 2)
        print("q=%d random: min (#nonzero + s' - (q+1)/2) = %d" % (q, minslack))
    print("Part 1 vectors tested:", tested)


# ---------------------------------------------------------------- Part 2
def structured(H, rng_count):
    """yield random and structured zero-sum integer vectors."""
    M = H.shape[0]
    q = M - 1
    for _ in range(rng_count):
        t = random.randint(0, 9)
        c = np.zeros(M, dtype=np.int64)
        if t == 0:  # random sparse ternary
            s = random.randint(2, M // 2) & ~1
            supp = random.sample(range(M), s)
            c[supp[: s // 2]] = 1
            c[supp[s // 2:]] = -1
        elif t == 1:  # random integer
            s = random.randint(2, M)
            supp = random.sample(range(M), s)
            vals = [random.randint(-3, 3) for _ in range(s - 1)]
            vals.append(-sum(vals))
            c[supp] = vals
        elif t == 2:  # column restricted to a balanced subset
            a = random.randint(1, q)
            h = H[:, a]
            P = [x for x in range(M) if h[x] == 1]
            N_ = [x for x in range(M) if h[x] == -1]
            k = random.randint(1, M // 2)
            c[random.sample(P, k)] = 1
            c[random.sample(N_, k)] = -1
        elif t == 3:  # column plus small perturbation
            a = random.randint(1, q)
            c = H[:, a].copy()
            for _ in range(random.randint(1, 4)):
                x, y = random.sample(range(M), 2)
                c[x] += 1
                c[y] -= 1
        elif t == 4:  # half-sum plus perturbation
            a, b = random.sample(range(1, M), 2)
            c = (H[:, a] + random.choice([1, -1]) * H[:, b]) // 2
            for _ in range(random.randint(1, 4)):
                x, y = random.sample(range(M), 2)
                c[x] += 1
                c[y] -= 1
        elif t == 5:  # combination of a few columns, halved if possible
            k = random.randint(1, 4)
            cols = random.sample(range(1, M), k)
            v = sum(random.choice([1, -1]) * H[:, a] for a in cols)
            for d in (4, 2, 1):
                if (v % d == 0).all():
                    c = v // d
                    break
        elif t == 6:  # half-sum restricted to a balanced subset
            a, b = random.sample(range(1, M), 2)
            g = (H[:, a] + H[:, b]) // 2
            P = [x for x in range(M) if g[x] == 1]
            N_ = [x for x in range(M) if g[x] == -1]
            k = random.randint(1, min(len(P), len(N_)))
            c[random.sample(P, k)] = 1
            c[random.sample(N_, k)] = -1
        elif t == 7:  # sum of two pairs / small multiples
            x, y, z, w = random.sample(range(M), 4)
            c[x] += 1; c[y] -= 1; c[z] += random.choice([1, 2]); c[w] -= c[z]
        elif t == 8:  # column times 2 plus pair
            a = random.randint(1, q)
            c = 2 * H[:, a]
            x, y = random.sample(range(M), 2)
            c[x] += 1; c[y] -= 1
        else:  # dense random ternary
            perm = np.random.permutation(M)
            c[perm[: M // 2]] = 1
            c[perm[M // 2: M]] = -1
            c[perm[M - random.randint(0, M // 2):]] = 0
            if c.sum() != 0:
                c[perm[0]] -= c.sum()
        if c.any() and c.sum() == 0:
            yield c


def part2():
    print("== Part 2: identities and inequalities used in the proof ==")
    total = 0
    for q in (11, 19, 23, 43, 47, 59, 67, 83):
        H = paley(q)
        M = q + 1
        for c in structured(H, 4000):
            T = transform(H, c)
            n = int(np.abs(c).sum())
            m2 = int((c * c).sum())
            Fc = int(np.abs(T).sum())
            total += 1
            # (D)
            check(Fc * n == M * m2 + int((np.abs(T) * (n - np.abs(T))).sum()), "D q=%d" % q)
            # (I)
            k = int((2 * np.abs(T) > n).sum())
            check(abs(Fc - k * n) <= 2 * (Fc - M), "I q=%d" % q)
            # (L)
            check(Fc >= n and Fc >= M, "L q=%d" % q)
            # (P1)
            for a in range(1, M):
                if (c == H[:, a]).all() or (c == -H[:, a]).all():
                    continue
                Ta = abs(int(T[a - 1]))
                check(Fc >= min(2 * Ta, 2 * M), "P1 q=%d" % q)
            # (P2) for the two largest labels
            idx = np.argsort(-np.abs(T))[:2]
            a, b = idx[0] + 1, idx[1] + 1
            sgn = 1 if T[idx[0]] > 0 else -1
            cc = sgn * c
            Ta = int(sgn * T[idx[0]])
            Tb = int(sgn * T[idx[1]])
            sigma = 1 if Tb >= 0 else -1
            g2 = H[:, a] + sigma * H[:, b]
            if not (cc * 2 == g2).all():
                phi = lambda t: abs(t) - abs(t - M // 2)
                check(Fc >= M + phi(Ta) + phi(sigma * Tb), "P2 q=%d" % q)
            # (S)
            sp = int((c[:q] != 0).sum())
            if 4 <= n <= (q - 1) // 2:
                fn = (2 * (n - 2) / n) * (M / 2 - n - M / n)
                check(Fc - M >= fn - 1e-9, "S q=%d n=%d F-M=%d f=%.3f" % (q, n, Fc - M, fn))
    print("Part 2 vectors tested:", total)


# ---------------------------------------------------------------- Part 3
def hill(H, c, steps=60):
    """greedy descent of F over moves c -> c +- (e_x - e_y), staying non-exceptional and nonzero."""
    M = H.shape[0]
    Hn = H[:, 1:]
    T = transform(H, c)
    best = int(np.abs(T).sum())
    D = (Hn[:, None, :] - Hn[None, :, :]).reshape(M * M, M - 1)  # row x*M+y: H_x - H_y
    for _ in range(steps):
        Fs = np.abs(T[None, :] + D).sum(axis=1)
        order = np.argsort(Fs)
        moved = False
        for idx in order[:50]:
            if Fs[idx] >= best:
                break
            x, y = divmod(int(idx), M)
            if x == y:
                continue
            c2 = c.copy(); c2[x] += 1; c2[y] -= 1
            if not c2.any() or exceptional(H, c2):
                continue
            c, T, best = c2, T + D[idx], int(Fs[idx])
            moved = True
            break
        if not moved:
            break
    return c, best


def part3():
    print("== Part 3: search for near-flat non-exceptional characters ==")
    print("   (Lemma S claims min (F-M)/M >= 1/16 = 0.0625 among non-exceptional c, for q >= 43)")
    for q in (7, 11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83, 103, 107):
        H = paley(q)
        M = q + 1
        best = (10 ** 9, None)
        cnt = 0
        # exhaustive ternary of support 4 and 6 (support 6 only for small q), and {+-1,+-2} support 3
        maxs = 6 if q <= 23 else 4
        for s in range(3, maxs + 1):
            for supp in combinations(range(M), s):
                vals_list = (product([1, -1], repeat=s) if s % 2 == 0 else
                             product([1, -1, 2, -2], repeat=s))
                for vals in vals_list:
                    if sum(vals) != 0 or vals[0] < 0:
                        continue
                    c = np.zeros(M, dtype=np.int64)
                    c[list(supp)] = vals
                    cnt += 1
                    f = F(H, c)
                    if f - M < best[0] and not exceptional(H, c):
                        best = (f - M, c.copy())
            if cnt > 3_000_000:
                break
        ex_best = best[0]
        # structured + hill climbing
        starts = list(structured(H, 600))
        for c in starts:
            if exceptional(H, c):
                continue
            f = F(H, c)
            if f - M < best[0]:
                best = (f - M, c.copy())
        for c in starts[:120]:
            if exceptional(H, c):
                continue
            c2, f2 = hill(H, c, steps=40)
            cnt += 1
            if f2 - M < best[0]:
                best = (f2 - M, c2.copy())
        c = best[1]
        n = int(np.abs(c).sum())
        print("q=%3d M=%3d: exhaustive-small min F-M = %4d; overall min F-M = %4d  ((F-M)/M = %.4f), "
              "at ||c||_1 = %d, support %d" % (q, M, ex_best, best[0], best[0] / M, n,
                                               int((c != 0).sum())))
        if q >= 43:
            check(best[0] >= M / 16, "Lemma S violated at q=%d" % q)


if __name__ == "__main__":
    part1()
    part2()
    part3()
    print("ALL CHECKS PASSED" if not FAIL else "FAILURES: %d" % len(FAIL))
