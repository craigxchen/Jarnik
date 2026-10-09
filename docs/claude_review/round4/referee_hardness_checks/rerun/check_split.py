#!/usr/bin/env python3
"""Exact checks for Section 4 of hardness.md (the split twin).

S(n,D,C) = #{d | n : D <= d <= D + C D^(3/2) n^(-1/2)}   (both d and n/d in square-root windows).

(H) Coppersmith lattice (Lemma 4.4): for given (n, P, m, t, X) build the basis
      g_j(y) = n^(m-j)(y+P)^j (0<=j<=m),  g_j(y) = y^(j-m)(y+P)^m (m<j<omega),
    reduce it with exact LLL, take the first vector h, and verify (i) the triangular determinant
    formula, (ii) sum |h_i| X^i < n^(beta m)  [the Howgrave-Graham condition], and
    (iii) h(x0) = 0 for every x0 in [0,X] with P+x0 a divisor of n of size >= n^beta.
(I) Vertex: S(n, sqrt n, C) <= 1 + C^2 (exact, n <= 2*10^5 all n).
(J) Data: max over D in [n^(1/3), n^(1/2)] of S(n,D,1) for all n <= N0, split by the position
    lambda = sqrt(n)/D of the record window.
"""
import math, random, sys
from fractions import Fraction
random.seed(11)
out = []
def log(s):
    print(s); out.append(s)

def poly_mul(a, b):
    r = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i+j] += x*y
    return r
def poly_pow(a, e):
    r = [1]
    for _ in range(e): r = poly_mul(r, a)
    return r
def poly_eval(p, x):
    s = 0
    for c in reversed(p): s = s*x + c
    return s

def lll(B, delta=Fraction(3, 4)):
    B = [list(map(int, v)) for v in B]
    n = len(B)
    def dot(u, v): return sum(a*b for a, b in zip(u, v))
    def gso(B):
        Bs = []; mu = [[Fraction(0)]*n for _ in range(n)]
        for i in range(n):
            v = [Fraction(x) for x in B[i]]
            for j in range(i):
                mu[i][j] = Fraction(dot(B[i], Bs[j]))/dot(Bs[j], Bs[j]) if dot(Bs[j], Bs[j]) else Fraction(0)
                v = [a - mu[i][j]*b for a, b in zip(v, Bs[j])]
            Bs.append(v)
        return Bs, mu
    Bs, mu = gso(B)
    k = 1
    while k < n:
        for j in range(k-1, -1, -1):
            q = round(mu[k][j])
            if q:
                B[k] = [a - q*b for a, b in zip(B[k], B[j])]
                Bs, mu = gso(B)
        if dot(Bs[k], Bs[k]) >= (delta - mu[k][k-1]**2)*dot(Bs[k-1], Bs[k-1]):
            k += 1
        else:
            B[k], B[k-1] = B[k-1], B[k]
            Bs, mu = gso(B)
            k = max(k-1, 1)
    return B

def coppersmith_poly(n, P, m, t, X):
    omega = m + 1 + t
    polys = []
    for j in range(m+1):
        polys.append([n**(m-j)*c for c in poly_pow([P, 1], j)])
    for j in range(m+1, omega):
        polys.append([0]*(j-m) + poly_pow([P, 1], m))
    basis = []
    for p in polys:
        v = [ (p[i] if i < len(p) else 0) * X**i for i in range(omega)]
        basis.append(v)
    det = 1
    for j, v in enumerate(basis):
        det *= v[j]
        assert all(v[i] == 0 for i in range(j+1, omega))
    assert det == n**(m*(m+1)//2) * X**(omega*(omega-1)//2)
    red = lll(basis)
    hv = min(red, key=lambda v: sum(x*x for x in v))
    h = [Fraction(hv[i], X**i) for i in range(omega)]
    assert all(c.denominator == 1 for c in h)
    return [int(c) for c in h], det, omega

def divisors(n):
    f = {}; x = n; p = 2
    while p*p <= x:
        while x % p == 0:
            f[p] = f.get(p, 0)+1; x //= p
        p += 1
    if x > 1: f[x] = f.get(x, 0)+1
    ds = [1]
    for p, e in f.items():
        ds = [d*p**k for d in ds for k in range(e+1)]
    return sorted(ds)


def divisors_of_lcm(ds_, n):
    # divisors of n lying near P: test candidates directly (window is tiny)
    P = ds_[0]
    return [d for d in range(P, P + 25) if n % d == 0]

# (H)
log("(H) Coppersmith lattice: exact LLL, Howgrave-Graham condition, vanishing at divisors")
cases = 0; vanish = 0
for trial in range(10):
    P = random.randint(10**12, 10**13)
    d1, d2 = sorted(random.sample(range(1, 21), 2))
    ds_ = [P, P + d1, P + d2]
    n = 1
    for d in ds_: n = n*d//math.gcd(n, d)
    m, t, X = 3, 3, 24
    h, det, omega = coppersmith_poly(n, P, m, t, X)
    lhs = sum(abs(c)*X**i for i, c in enumerate(h))
    cond = lhs < P**m          # every divisor d in [P, P+X] has d^m >= P^m > |h(d-P)|
    divs = [d for d in divisors_of_lcm(ds_, n) if P <= d <= P + X]
    for d in divs:
        assert poly_eval(h, d - P) % d**m == 0
    roots = [x for x in range(0, X+1) if poly_eval(h, x) == 0]
    if cond:
        assert all(poly_eval(h, d - P) == 0 for d in divs)
        assert len(roots) <= omega - 1
        vanish += len(divs)
    cases += 1
    log(f"  n~10^{math.log10(n):.1f}, P~10^{math.log10(P):.1f} (=n^{math.log(P)/math.log(n):.3f}), "
        f"omega={omega}, X={X}: HG condition {cond}; divisors of n in [P,P+X]: {len(divs)}; "
        f"integer roots of h in [0,X]: {roots}")
log(f"  {cases} cases; {vanish} divisor offsets certified as exact roots of one polynomial each")

# (I) vertex
log("(I) vertex: #{d | n : sqrt n <= d <= sqrt n + C n^(1/4)} <= 1 + C^2")
for C in [1, 2, 3]:
    mx = 0
    for n in range(2, 200001):
        r = math.isqrt(n)
        lo = r if r*r == n else r + 1
        hi_f = math.sqrt(n) + C*n**0.25
        cnt = sum(1 for d in range(lo, int(hi_f)+1) if n % d == 0)
        mx = max(mx, cnt)
    assert mx <= 1 + C*C
    log(f"  C={C}: max over n <= 2*10^5 is {mx} <= {1+C*C}")

# (J) data
N0 = int(sys.argv[1]) if len(sys.argv) > 1 else 400000
log(f"(J) max_D S(n,D,1), D in [n^(1/3), n^(1/2)], all n <= {N0}")
spf = list(range(N0+1))
for p in range(2, int(N0**0.5)+1):
    if spf[p] == p:
        for q in range(p*p, N0+1, p):
            if spf[q] == q: spf[q] = p
def divs_spf(n):
    ds = [1]
    while n > 1:
        p = spf[n]; e = 0
        while n % p == 0: n //= p; e += 1
        ds = [d*p**k for d in ds for k in range(e+1)]
    return sorted(ds)
best = {}  # bucket by lambda
records = []
for n in range(2, N0+1):
    ds = divs_spf(n)
    sq = math.sqrt(n); lo = n**(1/3)
    cand = [d for d in ds if lo <= d <= sq]
    for i, d in enumerate(cand):
        L = d**1.5/sq
        cnt = 0
        for e in cand[i:]:
            if e <= d + L: cnt += 1
            else: break
        if cnt >= 3:
            lam = sq/d
            b = 'vertex lambda<1.0001' if lam < 1.0001 else ('1<lambda<2' if lam < 2 else 'lambda>=2')
            if cnt > best.get(b, (0,))[0]:
                best[b] = (cnt, n, d, lam)
for b, v in sorted(best.items()):
    log(f"  {b:>22}: max count {v[0]} (first at n={v[1]}, D={v[2]}, lambda={v[3]:.5f})")

with open(__file__.replace('.py', '_output.txt'), 'w') as fh:
    fh.write('\n'.join(out) + '\n')
