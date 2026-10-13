"""Evidence/exact for Theorems A, B, C: Legendre patterns at q = 7 mod 8.
(1) Lemma 2.1: least prime q = 7 mod 8, q not | n, with (n/q) = -1, for sub-genuine n = p*Norm(zeta).
(2) Lemma 4.2: least such q with (n_c/q) = -1 for all of P pieces (exists by duality; size ~ 2^P expected).
(3) Theorem B: max parity agreement over many small zeta, per column, on Q_7 up to X."""
import random, math, sys
sys.path.insert(0, '.')
from tq import is_prime, primes_upto

def leg(a, q):
    t = pow(a % q, (q-1)//2, q)
    return 0 if t == 0 else (1 if t == 1 else -1)

rng = random.Random(7)
Q7 = [q for q in primes_upto(2*10**6) if q % 8 == 7]
def least_bad(ns):
    for q in Q7:
        if all(n % q != 0 and leg(n, q) == -1 for n in ns): return q
    return None

ps = []
while len(ps) < 200:
    p = rng.randrange(10**5, 10**7)
    if p % 4 == 1 and is_prime(p): ps.append(p)
norms = [a*a+b*b for a in range(1, 40) for b in range(0, 40) if a*a+b*b > 1]
# (1)
worst = 0
for p in ps:
    for _ in range(20):
        nz = rng.choice(norms)
        if nz % p == 0: continue
        n = p*nz
        s = math.isqrt(n); s2 = math.isqrt(n//2)
        if s*s == n or (n % 2 == 0 and 2*s2*s2 == n): continue
        q = least_bad([n]); worst = max(worst, q)
print(f"(1) Lemma 2.1: 4000 sub-genuine n = p*Norm(zeta), p in (1e5,1e7): max over n of least q = {worst}")
# (2)
for P in range(1, 13):
    vals = []
    for p in ps[:40]:
        ns = []
        while len(ns) < P:
            nz = rng.choice(norms)
            if nz % p: ns.append(p*nz)
        vals.append(least_bad(ns))
    vals = [v for v in vals if v]
    print(f"(2) P={P:2d} pieces: least common bad q over 40 columns: median {sorted(vals)[len(vals)//2]}, max {max(vals)}; log2(max) = {math.log2(max(vals)):.1f}")
# (3)
X = 10**6
Q7X = [q for q in Q7 if q <= X]; th = sum(math.log(q) for q in Q7X)
best = []
for p in ps[:30]:
    b = 0
    for nz in set(norms):
        if nz % p == 0: continue
        n = p*nz
        S = sum(math.log(q)*leg(n, q) for q in Q7X if n % q)
        b = max(b, S/th)
    best.append(b)
print(f"(3) X=1e6: over 30 columns and {len(set(norms))} norms each, max relative bias S_n/theta_7 = {max(best):.4f} (agreement fraction <= {(1+max(best))/2:.4f})")
