"""Lemma 3.1 of round6/sharp.md: the prime assignment q -> p'(q) and the pair demand of the
parity-respecting structured design.

Assignment.  p'(3) = p'(5) = 2, p'(7) = 3.  For k >= 3 let B_k = odd primes in [2^k, 2^(k+1)) and
T_k = primes in [2^(k-2), 2^(k-1)), listed increasingly; the i-th element of B_k (i = 0..n-1)
goes to the floor(i t/n)-th element of T_k (t = |T_k|).  Then p'(q) < 2^(k-1) <= (q-1)/2 <= m_q/2
and q < 8 p'(q).

Part 1 (exact, sieve up to 2^KMAX): multiplicities and Lambda'(p)/log p, where
Lambda'(p) = sum of log q over the fibre of p.
Part 2 (explicit bound for k >= 20): with Rosser-Schoenfeld  x/ln x < pi(x) (x >= 17) and
pi(x) < 1.25506 x/ln x (x > 1):  n/t <= R(k); we evaluate R(k) for every 20 <= k <= 2000 and
bound k > 2000 analytically (R(k) <= 1.51013/(0.186235 - 1/2000)).
Part 3 (exact): max over 1 <= d < M of the weighted pair demand
   D(d) = sum_{q < 4M, p'(q) | d} (1 + min(v_q(d), a*(q) - 1)) log q,
   a*(q) = min{a : p'(q) q^(a-1) >= M},
which bounds log m_xy (line 1) for every pair, whatever the parity colouring."""
import math, sys
from gres import primes_upto

KMAX = 24




def assignment(kmax):
    P = primes_upto(2 ** (kmax + 1))
    odd = [q for q in P if q > 2]
    pp = {3: 2, 5: 2, 7: 3}
    for k in range(3, kmax + 1):
        B = [q for q in odd if 2 ** k <= q < 2 ** (k + 1) and q > 7]
        T = [p for p in P if 2 ** (k - 2) <= p < 2 ** (k - 1)]
        n, t = len(B), len(T)
        assert t >= 1
        for i, q in enumerate(B):
            pp[q] = T[(i * t) // n]
    return pp


if __name__ == "__main__":
    KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 24
    P = primes_upto(2 ** (KMAX + 1))
    odd = [q for q in P if q > 2]
    pp = assignment(KMAX)
    ok = True
    mult = {}
    lam = {}
    for q, p in pp.items():
        m = q - 1 if q % 4 == 1 else q + 1
        ok &= (2 * p <= m) and (q < 8 * p) and p != q
        mult[p] = mult.get(p, 0) + 1
        lam[p] = lam.get(p, 0.0) + math.log(q)
    maxmult = max(mult.values())
    worst = max(lam[p] / math.log(p) for p in lam)
    print(f"Part 1: odd primes q < 2^{KMAX+1}: {len(pp)};  max multiplicity {maxmult};  max Lambda'(p)/log p = {worst:.3f}")
    by_k = {}
    for p, c in mult.items():
        k = p.bit_length() + 1
        by_k[k] = max(by_k.get(k, 0), c)
    print("   max multiplicity by block k:", [by_k.get(k) for k in range(3, KMAX + 1)])


    # Part 2: n <= pi(2^(k+1)) - pi(2^k) + 1 and t >= pi(2^(k-1)) - pi(2^(k-2)) - 1, with the
    # Rosser-Schoenfeld bounds; everything divided by 2^k / ln 2.  Every k in [20, 2000] is
    # evaluated; for k > 2000: num <= 1.51013/k and den >= (0.186235 - 1/k - 1e-500)/k.
    Rmax = 0.0
    for k in range(20, 2001):
        eps = math.log(2) * 2.0 ** (-k)
        num = 1.25506 * 2 / (k + 1) - 1 / k + eps
        den = 0.5 / (k - 1) - 1.25506 * 0.25 / (k - 2) - eps
        Rmax = max(Rmax, num / den)
    tail = 1.51013 / (0.186235 - 1.0 / 2000)
    print(f"Part 2: max_(20<=k<=2000) R(k) = {Rmax:.4f};  bound for k > 2000: {tail:.4f};  so multiplicity <= {math.ceil(max(Rmax, tail))} for k >= 20")
    ok &= max(Rmax, tail) < 9


    # Part 3
    def vq(d, q):
        v = 0
        while d % q == 0:
            d //= q
            v += 1
        return v


    for M in (64, 256, 1024, 4096, 16384):
        qs = [q for q in odd if q < 4 * M]
        astar = {}
        for q in qs:
            a = 1
            while pp[q] * q ** (a - 1) < M:
                a += 1
            astar[q] = a
        worstD, argd = 0.0, None
        for d in range(1, M):
            D = 0.0
            for q in qs:
                if d % pp[q] == 0:
                    D += (1 + min(vq(d, q), astar[q] - 1)) * math.log(q)
            if D > worstD:
                worstD, argd = D, d
        print(f"Part 3: M = {M:6d}: max pair demand {worstD:8.2f} at d = {argd:6d};  = {worstD / math.log(M):.3f} log M")
    print("ALL DESIGN CHECKS PASSED" if ok else "FAILURE")
