# C6 (evidence): parity bias S_n / theta_7(X) of sub-genuine n = p Norm(zeta) over q = 7 mod 8, q <= X = 10^6,
# for 30 random column primes and all distinct norms a^2+b^2 (1 <= a < 40, 0 <= b < 40). Multiplicativity of
# the Legendre symbol: vectors (l/q) for each prime l are computed once (exact integers), then multiplied.
import numpy as np, random, math
from glib import *
X = 10 ** 6
qs = np.array([q for q in sieve(X) if q % 8 == 7], dtype=np.int64)
logs = np.log(qs.astype(float)); th7 = logs.sum()
cache = {}
def leg_vec(l):
    if l not in cache:
        if l == 2:
            cache[l] = np.ones(len(qs), dtype=np.int8)   # (2/q) = 1 for q = 7 mod 8
        else:
            cache[l] = np.array([0 if q == l else (1 if pow(l, (q - 1) // 2, q) == 1 else -1) for q in qs.tolist()], dtype=np.int8)
    return cache[l]
def vec(n):
    v = np.ones(len(qs), dtype=np.int8)
    for l, e in factor(n).items():
        if e % 2: v = v * leg_vec(l)
        else:
            v = v * (leg_vec(l) != 0)  # q | l: symbol 0 (excluded from S_n)
    return v
rng = random.Random(3)
cols = []
while len(cols) < 30:
    p = rng.randrange(10 ** 5, 10 ** 7)
    if p % 4 == 1 and is_prime(p): cols.append(p)
norms = sorted(set(a * a + b * b for a in range(1, 40) for b in range(0, 40)))
worst = 0.0
for p in cols:
    vp = leg_vec(p)
    for nz in norms:
        if nz % p == 0: continue
        S = float((vp * vec(nz) * logs).sum())
        worst = max(worst, S / th7)
print(f"{len(cols)} columns x {len(norms)} norms: max S_n/theta_7(10^6) = {worst:.4f} "
      f"(agreement in parity at most {(1 + worst) / 2:.4f} of the Q_7 log-mass)")
