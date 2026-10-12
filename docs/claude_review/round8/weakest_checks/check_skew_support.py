"""Skew double core (round5/exact.md Thm 4.1): support lemma for characters.

K = I + S, skew Paley of order M = q+1 (q = 3 mod 4 prime), rows/labels {inf} u F_q,
S(inf,a) = -1, S(a,inf) = 1, S(a,b) = chi(b-a).  Cuts (columns of the profile):
S_a = {x : K(x,a) = -1},  S'_a = S_a xor {inf, a},  a in F_q   (2q cuts).
For a zero-sum lambda the character value on a cut T is lambda(T) = sum_{x in T} lambda_x.

Claim (Lemma 3.2 of weakest.md): every nonzero zero-sum rational lambda has
#{T : lambda(T) != 0} >= (q+1)/2.

Checks: K is skew-Hadamard; cuts balanced; every pair of rows separated by exactly M cuts;
then the support lemma (exhaustive over small supports, random dense vectors, and the
structured vectors K e_a +- K e_b that come closest), reporting the minimum found.
"""
import itertools
import random


def legendre(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1


def skew_core(q):
    idx = ['inf'] + list(range(q))
    M = q + 1
    S = [[0] * M for _ in range(M)]
    for j in range(1, M):
        S[0][j] = -1
        S[j][0] = 1
    for i in range(1, M):
        for j in range(1, M):
            S[i][j] = legendre((j - 1) - (i - 1), q)
    K = [[(1 if i == j else 0) + S[i][j] for j in range(M)] for i in range(M)]
    return K


def cuts(K):
    M = len(K)
    T = []
    for a in range(1, M):
        Sa = frozenset(x for x in range(M) if K[x][a] == -1)
        T.append(Sa)
        T.append(Sa ^ frozenset({0, a}))
    return T


def check_structure(q):
    K = skew_core(q)
    M = q + 1
    for i in range(M):
        for j in range(M):
            d = sum(K[i][t] * K[j][t] for t in range(M))
            assert d == (M if i == j else 0)
            assert K[i][j] + K[j][i] == (2 if i == j else 0)
    T = cuts(K)
    assert len(set(T)) == 2 * q
    assert all(len(c) == M // 2 for c in T)
    for x in range(M):
        for y in range(x + 1, M):
            sep = sum(1 for c in T if (x in c) != (y in c))
            assert sep == M, (x, y, sep)
    return K, T


def nsupp(lam, T):
    return sum(1 for c in T if sum(lam[x] for x in c) != 0)


def run(q, max_supp, rand_trials, seed=1):
    K, T = check_structure(q)
    M = q + 1
    bound = (q + 1) // 2
    best = (10 ** 9, None)
    count = 0
    # exhaustive: zero-sum integer vectors with entries in {-2..2}, support <= max_supp
    for s in range(2, max_supp + 1):
        for supp in itertools.combinations(range(M), s):
            for vals in itertools.product((-2, -1, 1, 2), repeat=s - 1):
                last = -sum(vals)
                if last == 0 or abs(last) > 2:
                    continue
                lam = [0] * M
                for x, v in zip(supp, vals + (last,)):
                    lam[x] = v
                count += 1
                n = nsupp(lam, T)
                if n < best[0]:
                    best = (n, tuple(lam))
    # structured: lambda = K(e_a + s e_b) (columns a,b >= 1 are balanced, so zero-sum)
    for a in range(1, M):
        for b in range(1, M):
            if a == b:
                continue
            for sgn in (1, -1):
                lam = [K[x][a] + sgn * K[x][b] for x in range(M)]
                assert sum(lam) == 0
                n = nsupp(lam, T)
                count += 1
                if n < best[0]:
                    best = (n, ('K(e_a%se_b)' % ('+' if sgn > 0 else '-'), a, b))
    rng = random.Random(seed)
    for _ in range(rand_trials):
        s = rng.randint(2, M)
        supp = rng.sample(range(M), s)
        lam = [0] * M
        for x in supp[:-1]:
            lam[x] = rng.choice((-3, -2, -1, 1, 2, 3))
        lam[supp[-1]] = -sum(lam)
        if all(v == 0 for v in lam):
            continue
        n = nsupp(lam, T)
        count += 1
        if n < best[0]:
            best = (n, tuple(lam))
    ok = best[0] >= bound
    print(f"q={q:>3} M={M:>3} cuts={2*q:>3}: vectors tested={count:>8}, min #nonzero cuts={best[0]:>3}, "
          f"bound (q+1)/2={bound:>3}, min/(2q)={best[0]/(2*q):.4f}  {'PASS' if ok else 'FAIL'}  witness={best[1] if not isinstance(best[1], tuple) or len(best[1]) < 6 else 'vector'}")
    return ok


if __name__ == '__main__':
    allok = True
    for q, ms, rt in ((7, 6, 20000), (11, 5, 20000), (19, 4, 20000), (23, 3, 20000), (31, 3, 10000), (43, 2, 5000)):
        allok &= run(q, ms, rt)
    print("ALL PASS" if allok else "SOME FAIL")


def union_support_pairs(q):
    """Second part (Remark 3.6.3): union support of lambda_1 = K(e_a - e_b), lambda_2 = K(e_a - e_c)."""
    K = skew_core(q)
    T = cuts(K)
    M = q + 1
    a = 1
    vals = []
    for b in range(1, M):
        for c in range(b + 1, M):
            if a in (b, c):
                continue
            l1 = [K[x][a] - K[x][b] for x in range(M)]
            l2 = [K[x][a] - K[x][c] for x in range(M)]
            u = sum(1 for t in T if sum(l1[x] for x in t) != 0 or sum(l2[x] for x in t) != 0)
            vals.append(u)
    print(f"q={q:>3}: union support of K(e_a-e_b), K(e_a-e_c): min={min(vals)}, mean={sum(vals)/len(vals):.1f}, "
          f"max={max(vals)}; min/(2q)={min(vals)/(2*q):.4f}, mean/(2q)={sum(vals)/len(vals)/(2*q):.4f} (heuristic 3/8)")


if __name__ == '__main__':
    for q in (43, 67, 103, 151):
        union_support_pairs(q)


def identities(q, trials=2000, seed=3):
    """Third part (Lemma 3.4 steps 2-3): lambda(S_a) = -(K^T lambda)_a/2, lambda(S'_a) = (K lambda~)_a/2,
    (K^T lambda)_inf = (K lambda~)_inf = 0, and both block sums of squares equal M ||lambda||^2/4."""
    K = skew_core(q)
    M = q + 1
    rng = random.Random(seed)
    for _ in range(trials):
        lam = [rng.randint(-5, 5) for _ in range(M - 1)]
        lam.append(-sum(lam))
        rng.shuffle(lam)
        lt = lam[:]
        lt[0] = -lt[0]
        KT = [sum(K[x][a] * lam[x] for x in range(M)) for a in range(M)]
        Kt = [sum(K[a][x] * lt[x] for x in range(M)) for a in range(M)]
        assert KT[0] == 0 and Kt[0] == 0
        b1 = b2 = 0
        for a in range(1, M):
            Sa = [x for x in range(M) if K[x][a] == -1]
            Spa = set(Sa) ^ {0, a}
            v1 = sum(lam[x] for x in Sa)
            v2 = sum(lam[x] for x in Spa)
            assert 2 * v1 == -KT[a] and 2 * v2 == Kt[a]
            b1 += v1 * v1
            b2 += v2 * v2
        n2 = sum(x * x for x in lam)
        assert 4 * b1 == M * n2 and 4 * b2 == M * n2
    print(f"q={q:>3}: identities of Lemma 3.4 steps 2-3 verified on {trials} random zero-sum vectors: PASS")


if __name__ == '__main__':
    for q in (7, 11, 19, 23, 43):
        identities(q)
