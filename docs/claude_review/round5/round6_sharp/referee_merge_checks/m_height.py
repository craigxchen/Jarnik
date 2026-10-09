"""Merge referee, check H: the height condition of sharp.md Prop. 6.1 / 6.5.

H1 (exact).  For actual Gaussian primes pi = a + bi (p = 1 mod 4, p < P1) and for products of two
     primes, compute D_s = Y, X-Y, X, X+Y of A_v and the odd-part collision moduli at ALL odd
     primes q (equivalently X = infinity, the largest possible moduli):
       per-unit   m_(v,s) = odd part of |D_s| (primes != p_j)  -- Prop. 6.1 says m_(v,s) <= nu_s^-1 |A_v|
       aggregated m_v     = prod_q q^(max_s a_(q,s))          -- sharp.md Verdict 5 / Prop 6.5 say <= sqrt2 |A_v|
     Also with X = 3000 (only q <= 3000 counted).
H2 (Monte Carlo of Prop. 6.5's random fibre, abstract model as in Lemma 6.3(2)).  r = b*M non-twin-like
     columns, each with independent fair Legendre bits nu_j(q) at every odd prime q <= X = 3M; label
     ell_j(q^a) uniform on the parity class nu_j(q) mod m_(q^a) (independent lifts).  i = g^(m/4).
     Unit s hit at level a iff ell = s*m/4 mod m_(q^a).  Report, in units of log^2 M/loglog M,
       tau_unit = 2 max_(j,s) log(nu_s m_(e_j,s))  (the least tau allowed by the per-unit condition)
       tau_aggr = 2 max_j log(m_(e_j)/sqrt2)      (what the aggregated condition would demand)
"""
import math, sys
import numpy as np

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

def gauss_prime(p):
    a = 1
    while True:
        b2 = p - a * a
        b = math.isqrt(b2)
        if b * b == b2: return a, b
        a += 1

def oddpart_upto(n, X, forbid=()):
    """product of q^v_q(n) over odd primes q <= X, q not in forbid."""
    n = abs(n); out = 1
    while n % 2 == 0: n //= 2
    q = 3
    rem = n
    while q * q <= rem and q <= X:
        while rem % q == 0:
            rem //= q
            if q not in forbid: out *= q
        q += 2
    if rem > 1 and rem <= X and rem not in forbid: out *= rem
    return out

def factor_odd(n, X, forbid=()):
    n = abs(n); f = {}
    while n % 2 == 0: n //= 2
    q = 3
    while q * q <= n:
        while n % q == 0:
            n //= q
            if q <= X and q not in forbid: f[q] = f.get(q, 0) + 1
        q += 2
    if n > 1 and n <= X and n not in forbid: f[n] = f.get(n, 0) + 1
    return f

def height_data(A, forbid, X):
    x, y = A
    Ds = [y, x - y, x, x + y]
    nu = [1.0, 2 ** -0.5, 1.0, 2 ** -0.5]
    absA = math.hypot(x, y)
    facs = [factor_odd(D, X, forbid) for D in Ds]
    per = [math.prod(q ** e for q, e in f.items()) for f in facs]
    allq = set().union(*[f.keys() for f in facs])
    agg = math.prod(q ** max(f.get(q, 0) for f in facs) for q in allq)
    viol_unit = any(per[s] > absA / nu[s] * (1 + 1e-12) for s in range(4))
    viol_aggr = agg > math.sqrt(2) * absA * (1 + 1e-12)
    return viol_unit, viol_aggr, per, agg, absA

def H1(P1=200000, X=10 ** 18):
    ps = [p for p in primes_upto(P1) if p % 4 == 1]
    vu = va = 0
    for p in ps:
        a, b = gauss_prime(p)
        u, g, per, agg, absA = height_data((a, b), {p}, X)
        vu += u; va += g
    print(f"H1 single primes p<{P1} (X={'inf' if X > 1e17 else X}): {len(ps)} primes; per-unit violations {vu}; aggregated violations {va}")
    a, b = gauss_prime(13)
    print(f"   example pi = {a}+{b}i: D_s = {[b, a-b, a, a+b]}, aggregated m = {height_data((a,b),{13},X)[3]} vs sqrt2*sqrt13 = {math.sqrt(26):.3f}")
    # two-prime monomials pi_1 pi_2 and pi_1 conj(pi_2)
    rng = np.random.default_rng(1)
    sel = rng.choice(len(ps), size=(2000, 2), replace=True)
    vu = va = tot = 0
    for i, j in sel:
        if i == j: continue
        p1, p2 = ps[i], ps[j]
        a1, b1 = gauss_prime(p1); a2, b2 = gauss_prime(p2)
        for sgn in (1, -1):
            x = a1 * a2 - sgn * b1 * b2; y = a1 * sgn * b2 + b1 * a2
            u, g, *_ = height_data((x, y), {p1, p2}, X)
            vu += u; va += g; tot += 1
    print(f"H1 two-prime monomials: {tot}; per-unit violations {vu}; aggregated violations {va}")

def H2(M, b=33, trials=1, seed=0):
    rng = np.random.default_rng(seed)
    X = 3 * M
    qs = [q for q in primes_upto(X) if q > 2]
    r = b * M
    lm_unit = np.zeros((r, 4)); lm_aggr = np.zeros(r)
    for q in qs:
        chi = 1 if q % 4 == 1 else -1
        Aq = int(math.floor(math.log(X) / math.log(q) + 1e-12))
        nu = rng.integers(0, 2, size=r)          # Legendre bit of p_j at q
        # level 1 labels: uniform on parity class mod m
        m = q - chi
        lab = nu + 2 * rng.integers(0, m // 2, size=r)
        hit_s = np.full(r, -1)
        for s in range(4):
            hit_s[lab == (s * m // 4) % m] = s
        lev = np.where(hit_s >= 0, 1, 0)
        # lifts: given a hit at level a-1, the level-a label is hit with prob 1/q (uniform lift)
        for a in range(2, Aq + 1):
            up = (lev == a - 1) & (rng.random(r) < 1.0 / q)
            lev[up] = a
        idx = np.nonzero(lev)[0]
        lm_unit[idx, hit_s[idx]] += lev[idx] * math.log(q)
        lm_aggr[idx] += lev[idx] * math.log(q)     # one unit per prime, so max_s = the hit unit
    nu_s = np.array([1, 2 ** -0.5, 1, 2 ** -0.5])
    tau_unit = 2 * (lm_unit + np.log(nu_s)).max()
    tau_aggr = 2 * (lm_aggr.max() - 0.5 * math.log(2))
    unit = math.log(M) ** 2 / math.log(math.log(M))
    return tau_unit / unit, tau_aggr / unit

if __name__ == "__main__":
    H1()
    H1(X=3000)
    for M in (1000, 3000, 10000):
        res = [H2(M, seed=s) for s in range(3)]
        print(f"H2 M={M}: tau_unit/(log^2M/loglogM) = {[f'{u:.2f}' for u, _ in res]}; tau_aggr = {[f'{g:.2f}' for _, g in res]}")
