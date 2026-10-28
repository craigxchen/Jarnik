"""Check 3 (exact integer counts; evidence for Thm P(c), Lemma 1.3 and Lemma 1.2):
(a) lattice count lemma: |#{ideals: N <= y, angle in I} - |I| y/2| <= 2.53 y^(1/2) + 1.58 on random sectors (|I| <= pi/2);
(b) polar boxes: |count under multiplicative graft - genuine count| versus the proved bound 16.6 (1+tau') x^(1/2) tau_(1/2);
(c) Psi_S(X) versus binom(2s+L, L) and X^kappa, kappa = log(4 e s)/log P0."""
import math, numpy as np
from common import ideals, grafts
from sympy import primerange
from math import comb

X = 2_000_000
x, y, n = ideals(X)
S = [p for p in primerange(101, 2000) if p % 4 == 1][:6]
rng = np.random.default_rng(3)
phi = rng.uniform(0, math.pi / 2, size=len(S))
th, mult, mini, smooth, spart = grafts(x, y, n, S, phi)
Q = math.pi / 2
def count(ang, lo, w, n1, n2):
    d = np.mod(ang - lo, Q)
    return int(((d <= w) & (n >= n1) & (n <= n2)).sum())
worst_a = 0.0
for _ in range(400):
    yy = int(rng.integers(100, X)); lo = rng.uniform(0, Q); w = rng.choice([1e-3, 1e-2, 0.1, 1.0, Q * 0.999])
    c = count(th, lo, w, 1, yy); err = abs(c - w * yy / 2)
    worst_a = max(worst_a, err / (2.53 * math.sqrt(yy) + 1.58))
print(f"(a) max |count-area| / (2.53 y^1/2 + 1.58) over 400 random sectors: {worst_a:.4f}  (<=1 required)")
tau_half = np.prod([(1 - p ** -0.5) ** -2 for p in S]) - 1
taup = np.prod([(1 + p ** -0.5) ** 2 for p in S]) - 1
worst_b = 0.0; maxdiff = 0
for _ in range(400):
    xx = int(rng.integers(10000, X)); yy = int(rng.integers(1, xx)); lo = rng.uniform(0, Q); w = rng.choice([1e-3, 1e-2, 0.1, 1.0])
    dg = count(th, lo, w, xx - yy, xx); dm = count(mult, lo, w, xx - yy, xx)
    maxdiff = max(maxdiff, abs(dg - dm))
    worst_b = max(worst_b, abs(dg - dm) / (16.6 * (1 + taup) * math.sqrt(xx) * tau_half))
print(f"(b) S={S}, tau_1/2={tau_half:.4f}, tau'={taup:.4f}; max |graft-genuine| = {maxdiff}; max ratio to proved bound = {worst_b:.4f}")
# (c)
for (Xc, P0, s) in [(10**8, 101, 6), (10**10, 1009, 10), (10**12, 10007, 20)]:
    Sc = [p for p in primerange(P0, 100 * P0) if p % 4 == 1][:s]
    # count S-smooth ideals: exponent vectors over 2s Gaussian primes with sum v log p <= log X
    logs = [math.log(p) for p in Sc for _ in (0, 1)]
    L = int(math.log(Xc) / math.log(P0))
    def cnt(i, rem):
        if i == len(logs): return 1
        t = 0; v = 0
        while v * logs[i] <= rem + 1e-12:
            t += cnt(i + 1, rem - v * logs[i]); v += 1
        return t
    psi = cnt(0, math.log(Xc))
    kappa = math.log(4 * math.e * s) / math.log(P0)
    print(f"(c) X={Xc:.0e} P0={P0} s={s}: Psi_S={psi}, binom(2s+L,L)={comb(2*s+L, L)}, X^kappa={Xc**kappa:.3e} (kappa={kappa:.3f})")
