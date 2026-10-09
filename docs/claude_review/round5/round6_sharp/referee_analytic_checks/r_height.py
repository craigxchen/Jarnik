"""Referee check of the height condition (sharp.md Prop. 6.1, §0 item 5, Prop. 6.5).

Prop. 6.1 (correct, per unit):      m_(v,s) <= nu_s^(-1) e^(w(v)/2)   for each s in Z/4.
§0 item 5 / Prop. 6.5 (aggregated): m_v := prod_q q^(max_s a_(q,s)(v)) <= sqrt2 e^(w(v)/2).

For an actual monomial A_v = X + iY (here single Gaussian primes pi over p, and products of two),
D_s = Y, X - Y, X, X + Y and m_(v,s) = odd prime-power part of D_s supported on q^a <= Xmod.
Part 1: count actual Gaussian primes (p < 10^6) and actual 2-prime monomials violating each form,
        for Xmod = 10^7 (all moduli) and Xmod = 3 p^(1/2) (the P3 regime p_j ~ M^2, X = 3M).
Part 2: abstract independent-fibre model of Prop. 6.5 (labels uniform in their parity class at
        each odd prime q <= Xmod; hit iff label in the 4-torsion; the unit s of each hit recorded),
        r = 33 M columns, random Legendre parities nu_j(q).  Reports, as multiples of
        log^2 M/loglog M, the max over columns of the aggregated log m_(e_j) and of the per-unit
        max_s log m_(e_j, s), at Xmod = 3M.
"""
import math, random
import numpy as np

rng = random.Random(3)
g = np.random.default_rng(3)


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


SP = primes_upto(4000)


def odd_part_bounded(D, Xmod):
    """prod of q^a' over odd primes q | D, with a' = max{a <= v_q(D) : q^a <= Xmod}."""
    D = abs(D)
    while D % 2 == 0 and D > 0:
        D //= 2
    m = 1
    for q in SP:
        if q == 2:
            continue
        if q * q > D:
            break
        if D % q == 0:
            a = 0
            while D % q == 0:
                D //= q
                a += 1
            t = 1
            for _ in range(a):
                if t * q <= Xmod:
                    t *= q
            m *= t
    assert D < SP[-1] ** 2 or D == 1
    if D > 1 and D <= Xmod:
        m *= D            # remaining odd prime
    return m


def gauss_prime(p):
    x = math.isqrt(p)
    for a in range(1, x + 1):
        b2 = p - a * a
        b = math.isqrt(b2)
        if b * b == b2:
            return a, b


def per_unit(X, Y, Xmod):
    Ds = [Y, X - Y, X, X + Y]
    return [odd_part_bounded(D, Xmod) for D in Ds]


nu_inv = [1.0, math.sqrt(2), 1.0, math.sqrt(2)]
out = []
P = [p for p in primes_upto(10 ** 6) if p % 4 == 1]
for label, XmodF in (("Xmod = 1e7", lambda p: 10 ** 7), ("Xmod = 3 sqrt(p)", lambda p: 3 * math.isqrt(p))):
    viol_unit = 0
    viol_agg = 0
    viol_prod = 0
    tot = 0
    for p in P:
        a, b = gauss_prime(p)
        Xm = XmodF(p)
        ms = per_unit(a, b, Xm)
        A = math.sqrt(p)
        if any(ms[s] > nu_inv[s] * A + 1e-9 for s in range(4)):
            viol_unit += 1
        agg = 1
        for v in ms:
            agg *= v                      # coprime across s (one unit per prime q)
        if agg > math.sqrt(2) * A:
            viol_agg += 1
        if agg > A ** 4 / 4 + 1e-9:
            viol_prod += 1
        tot += 1
    out.append(f"Part 1 ({label}): {tot} actual Gaussian primes p < 1e6: per-unit height violations = {viol_unit}; "
               f"aggregated 'm_v <= sqrt2 e^(w/2)' violations = {viol_agg} ({100*viol_agg/tot:.1f}%); "
               f"m_v <= e^(2w)/4 violations = {viol_prod}")
# example
a, b = gauss_prime(13)
out.append(f"   example p = 13, pi = {a}+{b}i: D_s = {[b, a-b, a, a+b]}, m_(v,s) = {per_unit(a, b, 10**7)}, "
           f"m_v = {np.prod(per_unit(a, b, 10**7))} > sqrt(2*13) = {math.sqrt(26):.2f}; per-unit bounds nu_s^-1 sqrt13 = "
           f"{[round(x*math.sqrt(13),2) for x in nu_inv]}")
# 2-prime monomials pi1 * pi2 and pi1 * conj(pi2)
viol_unit = viol_agg = tot = 0
for it in range(20000):
    p1, p2 = rng.sample(P[:5000], 2)
    a1, b1 = gauss_prime(p1)
    a2, b2 = gauss_prime(p2)
    if rng.random() < 0.5:
        b2 = -b2
    X = a1 * a2 - b1 * b2
    Y = a1 * b2 + a2 * b1
    ms = per_unit(X, Y, 10 ** 7)
    A = math.sqrt(p1 * p2)
    if any(ms[s] > nu_inv[s] * A + 1e-9 for s in range(4)):
        viol_unit += 1
    if np.prod(ms) > math.sqrt(2) * A:
        viol_agg += 1
    tot += 1
out.append(f"Part 1b: {tot} actual 2-prime monomials: per-unit violations = {viol_unit}; aggregated violations = {viol_agg} ({100*viol_agg/tot:.1f}%)")

# Part 2: abstract Prop 6.5 model with per-unit bookkeeping
for M in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6):
    r = 33 * M
    Xmod = 3 * M
    agg = np.zeros(r)
    unit = np.zeros((4, r))
    for q in primes_upto(Xmod):
        if q == 2:
            continue
        m = q - 1 if q % 4 == 1 else q + 1
        # nu_j(q) random parity; 4-torsion {0, m/4, m/2, 3m/4}; iota = m/4 * (unit), parity of m/4
        # hit probability: if m = 4 mod 8: 2 torsion elements of each parity -> 4/m, s in {0,2} (even) or {1,3} (odd)
        # if m = 0 mod 8: all torsion even -> 8/m for even-parity labels (s any of 4), 0 for odd
        nh = g.binomial(r, min(1.0, 8.0 / m))
        cols = g.choice(r, size=nh, replace=False) if nh else np.array([], dtype=np.int64)
        lq = math.log(q)
        if m % 8 == 4:
            # each column hit w.p. 4/m: thin the 8/m sample by 1/2
            keep = g.random(nh) < 0.5
            cols = cols[keep]
            par = g.integers(0, 2, size=len(cols))           # nu_j(q)
            s = par + 2 * g.integers(0, 2, size=len(cols))    # two units of that parity
        else:
            par = g.integers(0, 2, size=len(cols))
            cols = cols[par == 0]                               # odd labels never hit
            s = g.integers(0, 4, size=len(cols))
        agg[cols] += lq
        np.add.at(unit, (s, cols), lq)
    scale = math.log(M) ** 2 / math.log(math.log(M))
    out.append(f"Part 2: M = {M:8d}: max_j log m_(e_j) [aggregated] = {agg.max()/scale:.3f}, max_j max_s log m_(e_j,s) [per unit] = "
               f"{unit.max()/scale:.3f}  (units of log^2 M/loglog M; Prop 6.5 asymptote 1 = tau/2 threshold at A = 1)")
print("\n".join(out))
