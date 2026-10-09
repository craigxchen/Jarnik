"""Referee (analytic lens): which height condition is necessary? (sharp.md Prop. 6.1 vs the
aggregated form used in the Verdict item 5 and in the proof of Prop. 6.5.)

For an ACTUAL Gaussian prime pi = X + iY (p = X^2 + Y^2 = 1 mod 4) and v = e_1, with residue data
r = pi mod q^(A_q) at every odd prime power q^a <= Xmod:
   per unit (Prop. 6.1):   m_(v,s) = prod_q q^(a_(q,s)) <= nu_s^(-1) sqrt(p)    (necessary, proved)
   aggregated:             m_v     = prod_q q^(max_s a_(q,s)) <= sqrt(2) sqrt(p)  (used in Prop 6.5)
a_(q,s) = max{a <= A_q : q^a | D_s},  D_s = Y, X - Y, X, X + Y (s = 0, 1, 2, 3).
Also 2-prime monomials v = e_1 + e_2 and e_1 - e_2 (A_v = pi1 pi2 or pi1 conj(pi2)).
Exact integer arithmetic.
"""
import math, random

def sieve(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

def two_squares(p):
    # Cornacchia-free brute force via Hermite-Serret: find x^2 = -1 mod p, then Euclid
    for g in range(2, p):
        x = pow(g, (p - 1) // 4, p)
        if x * x % p == p - 1:
            break
    a, b = p, x
    while b * b > p:
        a, b = b, a % b
    return b, math.isqrt(p - b * b)

def part_le(D, Xmod):
    """prod over odd primes q of q^(min(v_q(D), A_q)) with q^(A_q) <= Xmod (only q <= Xmod)."""
    D = abs(D)
    while D % 2 == 0:
        D //= 2
    res = 1
    q = 3
    n = D
    fac = {}
    while q * q <= n:
        while n % q == 0:
            fac[q] = fac.get(q, 0) + 1
            n //= q
        q += 2
    if n > 1:
        fac[n] = fac.get(n, 0) + 1
    for q, e in fac.items():
        a = 0
        while a < e and q ** (a + 1) <= Xmod:
            a += 1
        res *= q ** a
    return res

def per_unit_and_agg(X, Y, Xmod):
    Ds = [Y, X - Y, X, X + Y]
    ms = [part_le(D, Xmod) for D in Ds]
    agg = 1
    # aggregated: per prime the max exponent over s (ms are products of disjoint-ish prime sets)
    primes = set()
    for D in Ds:
        n = abs(D)
        q = 3
        while n % 2 == 0: n //= 2
        while q * q <= n:
            if n % q == 0:
                primes.add(q)
                while n % q == 0: n //= q
            q += 2
        if n > 1: primes.add(n)
    for q in primes:
        a = 0
        for m in ms:
            e = 0
            while m % q == 0:
                m //= q; e += 1
            a = max(a, e)
        agg *= q ** a
    return ms, agg

nu_inv = [1.0, math.sqrt(2), 1.0, math.sqrt(2)]
P = [p for p in sieve(200000) if p % 4 == 1]
for Xmod in (10 ** 12, 3 * 1000):
    viol_unit = viol_agg = 0
    for p in P:
        X, Y = two_squares(p)
        ms, agg = per_unit_and_agg(X, Y, Xmod)
        if any(ms[s] > nu_inv[s] * math.sqrt(p) + 1e-9 for s in range(4)):
            viol_unit += 1
        if agg > math.sqrt(2 * p) + 1e-9:
            viol_agg += 1
    print(f"single primes p < 2e5 ({len(P)}), moduli <= {Xmod:g}: per-unit (Prop 6.1) violations = {viol_unit};"
          f" aggregated 'm_v <= sqrt(2p)' violations = {viol_agg} ({100.0 * viol_agg / len(P):.1f}%)")
X, Y = two_squares(13)
ms, agg = per_unit_and_agg(X, Y, 10 ** 12)
print(f"example p = 13, pi = {X}+{Y}i: D_s = {[Y, X - Y, X, X + Y]}, m_(v,s) = {ms}, aggregated m_v = {agg} > sqrt(26) = {math.sqrt(26):.2f}")

rng = random.Random(5)
viol_unit = viol_agg = 0
N = 4000
for _ in range(N):
    p1, p2 = rng.sample(P[100:], 2)
    x1, y1 = two_squares(p1); x2, y2 = two_squares(p2)
    for sgn in (1, -1):
        X = x1 * x2 - sgn * y1 * y2
        Y = x1 * y2 * sgn + y1 * x2
        ms, agg = per_unit_and_agg(X, Y, 10 ** 12)
        w2 = math.sqrt(p1 * p2)
        if any(ms[s] > nu_inv[s] * w2 + 1e-9 for s in range(4)):
            viol_unit += 1
        if agg > math.sqrt(2) * w2 + 1e-9:
            viol_agg += 1
print(f"2-prime monomials ({2 * N}): per-unit violations = {viol_unit}; aggregated violations = {viol_agg} ({100.0 * viol_agg / (2 * N):.1f}%)")
