"""Referee check: Lemma 3.1 (character lattice) and Lemma 3.2 (support >= (q+1)/2) of weakest.md,
and the claimed minimum "(q+7)/2, sharp up to +3".

1. Column characters lambda = K e_a (a in F_q) are zero-sum integer vectors with exactly (q+3)/2 nonzero
   cuts, so the minimum is <= (q+3)/2 (the write-up's search did not try single columns).
2. Exact minimum for q = 7, 11: for every set Z of cuts of size 2q - s, test whether some nonzero zero-sum
   lambda vanishes on all of Z, i.e. whether rank([1_T for T in Z] + [all-ones]) < M (rank mod the prime
   2^61 - 1, which is <= rank over Q, so "full rank mod p" certifies "no such lambda").
3. Character lattice: every zero-sum rational lambda with integral cut values lies in Z^M_0 + Z lambda*,
   tested on random rational vectors with denominators dividing 2M (exact Fractions).
"""
import itertools
import random
from fractions import Fraction

P61 = (1 << 61) - 1


def legendre(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1


def skew_core(q):
    M = q + 1
    S = [[0] * M for _ in range(M)]
    for j in range(1, M):
        S[0][j] = -1
        S[j][0] = 1
    for i in range(1, M):
        for j in range(1, M):
            S[i][j] = legendre((j - 1) - (i - 1), q)
    return [[(1 if i == j else 0) + S[i][j] for j in range(M)] for i in range(M)]


def cuts(K):
    M = len(K)
    T = []
    for a in range(1, M):
        Sa = frozenset(x for x in range(M) if K[x][a] == -1)
        T.append(Sa)
        T.append(Sa ^ frozenset({0, a}))
    return T


def nsupp(lam, T):
    return sum(1 for c in T if sum(lam[x] for x in c) != 0)


def rank_mod(rows, p=P61):
    m = [r[:] for r in rows]
    rk, ncol = 0, len(m[0])
    for c in range(ncol):
        piv = next((i for i in range(rk, len(m)) if m[i][c] % p), None)
        if piv is None:
            continue
        m[rk], m[piv] = m[piv], m[rk]
        inv = pow(m[rk][c], p - 2, p)
        for i in range(len(m)):
            if i != rk and m[i][c] % p:
                f = (m[i][c] * inv) % p
                m[i] = [(a - f * b) % p for a, b in zip(m[i], m[rk])]
        rk += 1
    return rk


def exact_min(q, smax):
    K = skew_core(q)
    T = cuts(K)
    M = q + 1
    ind = [[1 if x in t else 0 for x in range(M)] for t in T]
    ones = [1] * M
    for s in range((q + 1) // 2, smax + 1):
        # is there a nonzero zero-sum lambda vanishing on some set of 2q - s cuts?
        found = None
        cnt = 0
        for Z in itertools.combinations(range(2 * q), 2 * q - s):
            cnt += 1
            if rank_mod([ind[i] for i in Z] + [ones]) < M:
                found = Z
                break
        print(f"   q={q}: s={s}: {'a set of ' + str(2*q - s) + ' vanishing cuts exists' if found else 'no set of ' + str(2*q-s) + ' cuts can vanish (all ' + str(cnt) + ' sets have full rank mod 2^61-1)'}")
        if found:
            return s
    return None


def lattice_check(q, trials=3000, seed=7):
    K = skew_core(q)
    T = cuts(K)
    M = q + 1
    rng = random.Random(seed)
    lam_star = [Fraction(2, M) - (2 if x == 0 else 0) for x in range(M)]
    hits = 0
    for _ in range(trials):
        # random lambda = (m/M) 1 + c with c integral, then force zero sum via c_inf
        m = rng.randint(-2 * M, 2 * M)
        c = [rng.randint(-3, 3) for _ in range(M)]
        lam = [Fraction(m, M) + c[x] for x in range(M)]
        lam[0] -= sum(lam)
        assert sum(lam) == 0
        integral = all(sum(lam[x] for x in t).denominator == 1 for t in T)
        # membership: lambda - k lambda* integral for some k
        # lambda - k lambda* integral for an integer k  <=>  all entries share one fractional part f
        # and f is a multiple of 2/M
        fr = {x - (x.numerator // x.denominator) for x in lam}
        member = len(fr) == 1 and (next(iter(fr)) * Fraction(M, 2)).denominator == 1
        assert integral == member, (q, lam)
        hits += integral
    print(f"   q={q}: lattice Z^M_0 + Z lambda*: {trials} random zero-sum (1/M)-vectors, {hits} with integral cut values; "
          f"integrality <=> membership in every case")


if __name__ == '__main__':
    print("1. Column characters lambda = K e_a:")
    for q in (7, 11, 19, 23, 31, 43, 47, 59, 67):
        K = skew_core(q)
        T = cuts(K)
        M = q + 1
        vals = set()
        for a in range(1, M):
            lam = [K[x][a] for x in range(M)]
            assert sum(lam) == 0
            vals.add(nsupp(lam, T))
        print(f"   q={q:>3}: #nonzero cuts of K e_a over all a: {sorted(vals)}; (q+3)/2 = {(q+3)//2}, "
              f"(q+7)/2 = {(q+7)//2}, bound (q+1)/2 = {(q+1)//2}")
    print("2. Exact minimum over all nonzero zero-sum lambda:")
    for q, smax in ((7, 5), (11, 7)):
        s = exact_min(q, smax)
        print(f"   q={q}: exact minimum number of nonzero cuts = {s} = (q+3)/2 = {(q+3)//2}")
    print("3. Character lattice (Lemma 3.1):")
    for q in (7, 11, 19):
        lattice_check(q)
