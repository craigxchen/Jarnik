"""Referee checks of sharp.md: Remark 2.3(2), Prop. 7.1 (Siegel), §7.3 count, Prop. 6.5 (per-unit).

(a) Remark 2.3(2): for actual Gaussian primes pi over p = 1 mod 4 and odd Q = q^a <= 300 with 8 | m,
    ell(pi) = log_g(pi/conj pi) mod 4 is unchanged under pi -> u pi (units) and rational rescaling,
    and both values {0,2} occur among p with (p/q) = +1 (so not a function of the Legendre symbol).
(b) Prop. 7.1 (Siegel): for M = 5..8 and random integer labels |lambda| <= B, a nonzero c with
    sum c = 0, <c, lambda> = 0 and ||c||_inf <= floor((M max(1,B))^(2/(M-2))) exists (exhaustive search).
(c) §7.3: exponential growth rate of #{c in Z^M : ||c||_1 <= 4M} (claimed ~ e^(5.3 M)).
(d) Prop. 6.5: per-unit height condition.  In the random-fibre model at q = +-3 mod 8 a column's
    hit at q has unit s of parity nu_j(q); Prop 6.1 constrains m_(v,s) for EACH s, while Prop 6.5
    bounds m_(e_j) = prod_q q^(max_s a).  We simulate (M up to 10^6) the max over columns of
    log m_(e_j) and of max_s log m_(e_j,s), in units of log^2 M/loglog M."""
import itertools, math, random
import numpy as np

rng = random.Random(99)


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def gm(u, v, Q):
    return ((u[0] * v[0] - u[1] * v[1]) % Q, (u[0] * v[1] + u[1] * v[0]) % Q)


def gp(u, e, Q):
    r_ = (1, 0)
    b_ = (u[0] % Q, u[1] % Q)
    while e:
        if e & 1:
            r_ = gm(r_, b_, Q)
        b_ = gm(b_, b_, Q)
        e >>= 1
    return r_


def mT(q, a):
    return (q - (1 if q % 4 == 1 else -1)) * q ** (a - 1)


# (a)
print("(a) quartic part of ell for actual Gaussian primes, 8 | m:")
oka = True
for q in [x for x in range(3, 300) if is_prime(x)]:
    a = 1
    while q ** (a + 1) <= 300:
        a += 1
    Q, m = q ** a, mT(q, a)
    if m % 8:
        continue
    # generator
    fs = [p for p in range(2, m + 1) if m % p == 0 and is_prime(p)]
    g = next((x, y) for x in range(Q) for y in range(Q)
             if (x * x + y * y) % Q == 1 and all(gp((x, y), m // p, Q) != (1, 0) for p in fs))
    tab, t = {}, (1, 0)
    for k in range(m):
        tab[t] = k
        t = gm(t, g, Q)
    vals = {1: set(), -1: set()}
    for p in [p for p in range(5, 3000) if p % 4 == 1 and is_prime(p) and p % q]:
        x = next(x for x in range(1, p) if math.isqrt(p - x * x) ** 2 == p - x * x)
        y = math.isqrt(p - x * x)
        ells = set()
        for u in [(1, 0), (0, 1), (-1, 0), (0, -1)]:
            for s in [1, 2, 3, 5]:
                if s % q == 0:
                    continue
                pi = gm(u, (s * x, s * y), Q)
                cpi = (pi[0], (-pi[1]) % Q)
                n = pow((cpi[0] ** 2 + cpi[1] ** 2) % Q, -1, Q)
                icpi = (cpi[0] * n % Q, (-cpi[1]) * n % Q)
                ells.add(tab[gm(pi, icpi, Q)] % 4)
        oka &= len(ells) == 1
        L = 1 if pow(p, (q - 1) // 2, q) == 1 else -1
        vals[L] |= ells
    oka &= vals[1] == {0, 2} and vals[-1] == {1, 3}
print("   invariant under units/rational rescaling and both quartic classes occur:", oka)

# (b)
print("(b) Siegel's lemma in Prop 7.1:")
okb = True
for M in range(5, 9):
    for B in ([1, 3, 10, 30, 100] if M <= 6 else [1, 3, 10, 30]):
        bound = int(math.floor((M * max(1, B)) ** (2 / (M - 2)) + 1e-12))
        for _ in range(30 if M <= 6 else 8):
            lam = [rng.randint(-B, B) for _ in range(M)]
            found = False
            for c in itertools.product(range(-bound, bound + 1), repeat=M - 1):
                cl = -sum(c)
                if abs(cl) > bound or (not any(c) and cl == 0):
                    continue
                cc = list(c) + [cl]
                if sum(ci * li for ci, li in zip(cc, lam)) == 0:
                    found = True
                    break
            okb &= found
    print(f"   M = {M}: all trials found c within the Siegel bound: {okb}")

# (c)
print("(c) #{c in Z^M : ||c||_1 <= 4M}:")
def logcount(M, L):
    # exact count via sum_k 2^k C(M,k) C(L,k)  (points of the l1 ball of radius L in Z^M)
    tot = sum(2 ** k * math.comb(M, k) * math.comb(L, k) for k in range(0, min(M, L) + 1))
    return math.log(tot)
for M in [50, 200, 800]:
    print(f"   M = {M}: (1/M) log count = {logcount(M, 4 * M) / M:.3f}")
# zero-sum constraint changes the rate by o(1)
best = max(b * math.log(2) + (-(b * math.log(b) + (1 - b) * math.log(1 - b)) if 0 < b < 1 else 0)
           + 4 * (-(b / 4 * math.log(b / 4) + (1 - b / 4) * math.log(1 - b / 4))) for b in np.linspace(0.01, 1, 991))
print(f"   entropy limit max_beta [beta ln2 + H(beta) + 4 H(beta/4)] = {best:.3f}  (claimed 5.3)")

# (d)
print("(d) Prop 6.5: max over columns of log m(e_j) (max over units) vs max_s log m(e_j, s)")
def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


nrng = np.random.default_rng(5)
for M in [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6]:
    P = primes_upto(M)
    Q0 = P[(P > 2) & ((P % 8 == 3) | (P % 8 == 5))]
    Q0 = Q0[Q0 >= M ** 0.5]
    ncol = 31 * M if M <= 10 ** 4 else 31 * 10 ** 4      # (b-2) M columns (capped)
    lq = np.log(Q0.astype(float))
    m = Q0 + np.where(Q0 % 4 == 3, 1, -1)
    best_any = best_unit = 0.0
    for start in range(0, ncol, 2000):
        nc = min(2000, ncol - start)
        # hit at q with prob 4/m; given a hit the unit is uniform on the 2 units of the column's
        # parity at q; parity nu_j(q) is +-1 with prob 1/2 (actual Legendre symbols of p_j vary)
        hit = nrng.random((nc, len(Q0))) < (4.0 / m)
        par = nrng.integers(0, 2, (nc, len(Q0)))
        which = nrng.integers(0, 2, (nc, len(Q0)))
        unit = 2 * which + par                         # s in {0,2} if par = 0, {1,3} if par = 1
        any_ = (hit * lq).sum(1)
        per = np.stack([((hit & (unit == s)) * lq).sum(1) for s in range(4)], 1).max(1)
        best_any = max(best_any, any_.max())
        best_unit = max(best_unit, per.max())
    sc = math.log(M) ** 2 / math.log(math.log(M))
    print(f"   M = {M:>8}: max log m(e_j) = {best_any/sc:.3f}, max_s log m(e_j,s) = {best_unit/sc:.3f}  [units of log^2 M/loglog M]")
print("ALL MISC CHECKS PASSED" if (oka and okb) else "FAILURE")
