"""Referee check (independent code): genuine labels obey the height bound per unit for v = e_j:
   prod_{q: q | D_s(pi_j)} q  <=  |D_s(pi_j)| <= sqrt2 * sqrt(p_j).
   Reproduces the prime set of the author's designs.py (seed 1, p_j ~ M^3) WITHOUT importing it,
   and computes the per-unit demand of e_j over q in [M^0.5, 3M] by direct divisibility."""
import sys, math, random
from sympy import isprime, primerange
M = int(sys.argv[1]) if len(sys.argv) > 1 else 50
r = M; Cexp = 3.0; seed = 1
random.seed(seed)
X = 3 * M
lo = int(M ** Cexp); P = set()
while len(P) < r:
    c = random.randrange(lo, 2 * lo)
    if c % 4 == 1 and c > X and isprime(c): P.add(c)
P = sorted(P)
def two_sq(p):
    for a in range(1, math.isqrt(p) + 1):
        b2 = p - a*a; b = math.isqrt(b2)
        if b*b == b2: return a, b
qs = [q for q in primerange(3, X + 1) if q >= M ** 0.5]
lM = math.log(M)
worst = 0; wj = None
for p in P:
    a, b = two_sq(p)
    D = [b, a - b, a, a + b]
    for s in range(4):
        dem = sum(math.log(q) for q in qs if D[s] % q == 0)
        bound = math.log(abs(D[s]))
        assert dem <= bound + 1e-9
        if dem > worst: worst, wj = dem, (p, a, b, s)
print(f"M={M}: max over j,s of genuine per-unit demand of e_j = {worst/lM:.3f} log M  (at {wj});"
      f"  a-priori cap log(sqrt(2*2M^3))/log M = {math.log(math.sqrt(2*2*M**3))/lM:.3f}")
