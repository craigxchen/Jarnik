"""Referee (analytic lens): pair demand of design D_A (sharp.md Lemma 3.1, Lemma 3.3(1)) and the
constant mu_M = 10 + 26/loglog M used in Theorem 5.1 step 4.

Own implementation of the dyadic assignment p'(q) (Lemma 3.1).
 (1) exact multiplicity of p' for q < 2^23, and Lemma 3.1(3): Lambda'(d) <= 9 log d + 18.72 omega(d)
     for all d < 2^20;
 (2) exact worst-case pair demand  D(d) = sum_{q<4M, p'(q)|d} min(a*(q), A_q, 1+v_q(d)) log q
     (parity ignored = worst case), max over d < M, against mu_M log M and 10 log M;
 (3) the chain written in sharp.md, (10 + 26/loglog d) log d, is NOT monotone in d and exceeds
     mu_M log M for small d; the correct route is omega(d) <= 1.3841 log M/loglog M for all d < M;
 (4) where the Theorem 5.1 pair union bound closes when the PROVED mu_M is used (b = 33, A <= 2):
     total <= (M^2/2) 4 (2 + pi/C) M^mu_M e^(-(b-4) tau/4 + 1/4), tau >= log(21 r^2 (log M)^2).
"""
import math
import numpy as np

LIM = 1 << 23
sieve = np.ones(LIM + 1, dtype=bool); sieve[:2] = False
for i in range(2, int(LIM ** 0.5) + 1):
    if sieve[i]:
        sieve[i * i::i] = False
primes = np.nonzero(sieve)[0]

pprime = {3: 2, 5: 2, 7: 3}
k = 3
while (1 << (k + 1)) <= LIM:
    qs = [int(q) for q in primes[(primes >= (1 << k)) & (primes < (1 << (k + 1)))] if q > 7]
    ts = [int(t) for t in primes[(primes >= (1 << (k - 2))) & (primes < (1 << (k - 1)))]]
    n, t = len(qs), len(ts)
    for i, q in enumerate(qs):
        pprime[q] = ts[(i * t) // n]
    k += 1
from collections import Counter
mult = Counter(pprime.values())
print(f"(1) q < 2^23: {len(pprime)} odd primes assigned; max multiplicity = {max(mult.values())}")
ok = True
for q, p in pprime.items():
    m = q - 1 if q % 4 == 1 else q + 1
    ok &= (2 * p <= m) and (q < 8 * p)
print(f"    2 p'(q) <= m_q and q < 8 p'(q) for all: {ok}")

# Lambda'(d) for d < 2^20 via preimage lists
D = 1 << 20
pre = {}
for q, p in pprime.items():
    pre.setdefault(p, []).append(q)
lam = np.zeros(D, dtype=float)
omega = np.zeros(D, dtype=np.int64)
for p in primes[primes < D]:
    p = int(p)
    omega[p::p] += 1
    s = sum(math.log(q) for q in pre.get(p, []))
    lam[p::p] += s
d = np.arange(D, dtype=float)
viol = (lam[2:] > 9 * np.log(d[2:]) + 18.72 * omega[2:] + 1e-9).sum()
print(f"    Lemma 3.1(3) violations for 2 <= d < 2^20: {int(viol)}")

# (2) exact worst-case pair demand
def vq(x, q):
    c = 0
    while x % q == 0:
        x //= q; c += 1
    return c

def demand_max(M, A=1):
    X = 3 * M ** A
    best, arg = 0.0, None
    small = [q for q in pprime if q < 4 * M]
    astar = {}
    for q in small:
        a = 1
        while pprime[q] * q ** (a - 1) < M:
            a += 1
        Aq = int(math.floor(math.log(X) / math.log(q) + 1e-12))
        while q ** (Aq + 1) <= X: Aq += 1
        while q ** Aq > X: Aq -= 1
        astar[q] = min(a, Aq)
    byp = {}
    for q in small:
        byp.setdefault(pprime[q], []).append(q)
    for dd in range(1, M):
        tot = 0.0
        x = dd
        for p, qs in byp.items():
            if dd % p == 0:
                for q in qs:
                    tot += min(astar[q], 1 + vq(dd, q)) * math.log(q)
        if tot > best:
            best, arg = tot, dd
    return best, arg

print("(2) exact worst-case pair demand max_d D(d) vs mu_M log M:")
for M in (64, 256, 1024, 4096, 16384, 65536):
    best, arg = demand_max(M)
    L = math.log(M)
    mu = 10 + 26 / math.log(L)
    print(f"    M = {M:6d}: max D = {best:7.2f} = {best / L:5.2f} log M (at d = {arg}); mu_M = {mu:5.2f}; bound holds: {best <= mu * L}")

# (3) non-monotone chain
print("(3) the d-wise form (10 + 26/loglog d) log d:")
for dd in (3, 4, 5, 8, 16, 30, 100, 1000):
    print(f"    d = {dd:5d}: {(10 + 26 / math.log(math.log(dd))) * math.log(dd):8.1f}", end="")
    print(f"   vs mu_M log M at M = 44: {(10 + 26 / math.log(math.log(44))) * math.log(44):6.1f}")
# correct route: omega(d) <= 1.3841 log M / loglog M for all 2 <= d < M (checked for M <= 2^20 at powers of 2)
okw = True
for e in range(5, 21):
    M = 1 << e
    okw &= omega[2:M].max() <= 1.3841 * math.log(M) / math.log(math.log(M))
print(f"    max_(d<M) omega(d) <= 1.3841 log M/loglog M for M = 2^5..2^20: {okw}")

# (4) where the pair union bound closes with the proved mu_M (b = 33, C = 1)
b, C = 33, 1.0
print("(4) Thm 5.1 pair term with PROVED mu_M (natural log of the bound), b = 33, C = 1, A <= 2:")
for L in (10, 100, 1e3, 1e4, 3e4, 4e4, 1e5):
    ll = math.log(L)
    mu = 10 + 26 / ll
    r_log = math.log(b) + L               # log r ~ log b + log M
    tau = math.log(21) + 2 * r_log + 2 * ll
    logtot = 2 * L - math.log(2) + math.log(4 * (2 + math.pi / C)) + mu * L - (b - 4) * tau / 4 + 0.25
    print(f"    log M = {L:8.0f}: mu_M = {mu:6.3f}; ln(total) = {logtot:12.1f}")
