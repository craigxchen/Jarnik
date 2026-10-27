"""Theorem P triple condition, all moduli charged (X >= sqrt(2) max |Z_v|): P0 = {v : 1 <= ||v||_1 <= K,
gcd(Re(Z_v^2), N) = 1} (Re Z^2 = D_1 D_3 is odd).  Additive triples {v, w, v+w} in P0 with ||v||_1 = 1 or 2 and
v + w in P0 (norm <= K).  Reports |P0|, |T| and [Z^r : <T>] (exact, Smith normal form via sympy for r <= 40, else
via a modular Hermite reduction).  Prime sets: the first r split primes, or rand / bal / last style sets.
Usage: python3 triple_allnorm.py K p1..pr   |   python3 triple_allnorm.py K first r"""
import sys, math, random
from sympy import Matrix, ZZ, primerange
from sympy.matrices.normalforms import smith_normal_form
from sllib import two_sq, Zv, VK
K = int(sys.argv[1]); import os; Ysm = int(os.environ.get("YSMOOTH", "0"))
if sys.argv[2] == 'first':
    r = int(sys.argv[3]); ps = [p for p in primerange(5, 10**6) if p % 4 == 1][:r]
else:
    ps = [int(x) for x in sys.argv[2:]]; r = len(ps)
pis = [two_sq(p) for p in ps]; Nrad = math.prod(ps)
from sympy import factorint
def smooth(x, Y):
    x = abs(x)
    for p in primerange(2, Y+1):
        while x % p == 0: x //= p
        if x == 1: return True
    return x == 1
def inP0(v):
    X, Y = Zv(pis, v); R2 = X*X - Y*Y
    return math.gcd(R2, Nrad) == 1 and (Ysm == 0 or smooth(R2, Ysm))
small = [v for v in VK(r, min(2, K)) if inP0(v)]
smallset = set(small)
T = set()
# triples (v, w, v+w) with ||v|| <= 2, ||w|| <= K-||v||... enumerate w over P0 of norm <= K-1 lazily: use small (norm<=2) only
for v in small:
    for w in small:
        u = tuple(a+b for a, b in zip(v, w))
        if any(u) and sum(map(abs, u)) <= K and inP0(u):
            T.update([v, w, u])
def index(vecs):
    # gcd of maximal minors via Hermite-like integer row reduction (exact)
    rows = [list(v) for v in vecs]
    piv_rows = []; col = 0
    A = rows
    det = 1
    for c in range(r):
        # gcd-reduce column c among remaining rows
        while True:
            nz = [row for row in A if row[c] != 0]
            if not nz: return 0
            mn = min(nz, key=lambda row: abs(row[c]))
            newA = []
            for row in A:
                if row is mn: continue
                if row[c] != 0:
                    f = row[c] // mn[c]
                    row = [a - f*b for a, b in zip(row, mn)]
                newA.append(row)
            if all(row[c] == 0 for row in newA):
                det *= abs(mn[c]); A = newA; break
            A = newA + [mn]
    return det
idx = index(sorted(T)) if T else 0
print(f"Ysmooth={Ysm} K={K} r={r} primes {ps[:6]}{'...' if r > 6 else ''}: |P0 ∩ V_2|={len(small)} |T|={len(T)} index[Z^r:<T>]={idx if idx else 'rank<r'}")
