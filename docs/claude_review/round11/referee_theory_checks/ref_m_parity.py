"""Referee check (exact big integers): the round-10 order m(n) and its parity.

m(n) = largest integer m with 8*binom(m+n-2, n-1) <= (2n+1)^(n-1)   (T = n).
The widened slab |<k,beta_p>| <= m/e_p of theory.md Sec 3.2 means |l_p(k)| <= m/2 with l_p integer,
so l_p takes at most 2*floor(m/2)+1 values on A_p: m values if m is odd (level-sum lemma kills UM_p).

Outputs: parity statistics of m(n) on [nlo, nhi], parity at skew-Hadamard (Paley q+1, q=3 mod 4 prime)
orders in the regime of the two barrier profiles, and the inequalities needed by the parity repair
(m_e = 2 floor(m/2) >= max(0.7 n(n-1), 4w)) and by the twin argument (w + 4T < m).
"""
import sys
from math import comb, lgamma, log, sqrt, floor

def m_exact(n):
    rhs = (2 * n + 1) ** (n - 1)
    # initial guess from lgamma
    target = (n - 1) * log(2 * n + 1) - log(8)
    def lg(m):
        return lgamma(m + n - 1) - lgamma(n) - lgamma(m)  # log binom(m+n-2, n-1)
    lo, hi = 1, 4 * n * n
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if lg(mid) <= target:
            lo = mid
        else:
            hi = mid
    m = lo
    # exact correction
    while 8 * comb(m + n - 2, n - 1) > rhs:
        m -= 1
    while 8 * comb(m + 1 + n - 2, n - 1) <= rhs:
        m += 1
    return m

def is_prime(q):
    if q < 2:
        return False
    i = 2
    while i * i <= q:
        if q % i == 0:
            return False
        i += 1
    return True

def w_of(n, r):
    return n * sqrt(8 * (n - 1) * log(4 * r))

if __name__ == '__main__':
    nlo, nhi = int(sys.argv[1]), int(sys.argv[2])
    odd = even = 0
    bad_fix = []
    ms = {}
    for n in range(nlo, nhi + 1):
        m = m_exact(n)
        ms[n] = m
        if m % 2:
            odd += 1
        else:
            even += 1
        me = 2 * (m // 2)
        # repair inequalities (r = n^2 as a generous column count; also r = 9n-8)
        for r in (n, 9 * n - 8, n * n):
            w = w_of(n, r)
            if (n - 1) >= 300 * log(4 * r):
                if not (me >= 0.7 * n * (n - 1) and me >= 4 * w and w + 4 * n < m):
                    bad_fix.append((n, r, m))
    print(f"n in [{nlo},{nhi}]: m odd for {odd} values, even for {even} values")
    print("first 30 (n, m, m mod 2):", [(n, ms[n], ms[n] % 2) for n in range(nlo, min(nhi, nlo + 29) + 1)])
    print("repair/twin inequality failures in regime (should be empty):", bad_fix[:10])
    # barrier profiles: Paley-I round-5 fake (r = 9M-8) and skew double core (r = 2(M-1)), M = q+1
    for name, rfun in (("round5 Paley fake  r=9M-8", lambda M: 9 * M - 8), ("skew double core r=2(M-1)", lambda M: 2 * (M - 1))):
        rows = []
        for M in range(nlo, nhi + 1):
            q = M - 1
            if q % 4 == 3 and is_prime(q) and (M - 1) >= 300 * log(4 * rfun(M)):
                rows.append((M, ms[M], ms[M] % 2))
        oddrows = [x for x in rows if x[2] == 1]
        print(f"{name}: Paley orders in regime within range: {len(rows)}; with m odd: {len(oddrows)}")
        print("   first 12:", rows[:12])
        print("   first 12 with m odd:", oddrows[:12])
