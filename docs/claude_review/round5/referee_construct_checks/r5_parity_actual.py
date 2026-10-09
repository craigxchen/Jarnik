"""Referee check R5: validate the prime-level parity law used in R4(c) on ACTUAL lattice points.

Claim: for actual z_x = eps_x prod_j pi_j^(a_xj) conj(pi_j)^(e_j - a_xj) on x^2+y^2 = N and odd prime q not dividing
N, the norm-one ratio z_x/z_0 mod q has class in K/K^2 (K = norm-one subgroup of (Z[i]/q)^*, cyclic of even
order q - chi(q)) equal to class(eps_x/eps_0) * prod_j chi_q(p_j)^(a_xj - a_0j), where class(i) = (2/q) and
class(-1) = +1.  Checked exactly for every pair of primitive points of several circles and every odd q <= 60
coprime to N (with eps_x computed from the factorisation)."""
import math, itertools

def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n**0.5) + 1))
def two_sq(p):
    for a in range(1, int(p**0.5) + 1):
        b = math.isqrt(p - a*a)
        if b*b == p - a*a: return (a, b)
def gmul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def gmulm(a, b, m): return ((a[0]*b[0] - a[1]*b[1]) % m, (a[0]*b[1] + a[1]*b[0]) % m)
def gconj(a): return (a[0], -a[1])
def gnorm(a): return a[0]*a[0] + a[1]*a[1]
def gpowm(a, k, m):
    r = (1, 0); a = (a[0] % m, a[1] % m)
    while k:
        if k & 1: r = gmulm(r, a, m)
        a = gmulm(a, a, m); k >>= 1
    return r
def ginvm(a, m):
    ni = pow(gnorm(a) % m, -1, m); c = gconj(a); return ((c[0]*ni) % m, (c[1]*ni) % m)
def gdiv(z, d):
    n = gnorm(d); w = gmul(z, gconj(d))
    return (w[0]//n, w[1]//n) if w[0] % n == 0 and w[1] % n == 0 else None
def legendre(a, q):
    a %= q; return 0 if a == 0 else (1 if pow(a, (q-1)//2, q) == 1 else -1)

def cls(y, q):
    """class of a norm-one residue y mod q in K/K^2: y^(|K|/2) = +-1."""
    order = q + 1 if q % 4 == 3 else q - 1
    t = gpowm(y, order // 2, q)
    assert t in ((1, 0), (q - 1, 0)), t
    return 1 if t == (1, 0) else -1

UNIT = {(1, 0): 0, (0, 1): 1, (-1, 0): 2, (0, -1): 3}
total = 0
for primes in [(5, 13, 17, 29), (5, 13, 17, 37, 41), (13, 17, 29, 53, 61)]:
    N = 1
    for p in primes: N *= p
    pis = [two_sq(p) for p in primes]
    pts = []
    for x in range(-math.isqrt(N), math.isqrt(N) + 1):
        y = math.isqrt(N - x*x)
        if y*y == N - x*x:
            pts += [(x, y)] + ([(x, -y)] if y else [])
    data = []
    for z in pts:
        a = []; w = z; ok = True
        for pi in pis:
            w1 = gdiv(w, pi)
            if w1 is not None and gdiv(w, gconj(pi)) is not None: ok = False; break
            if w1 is not None: a.append(1); w = w1
            else: a.append(0); w = gdiv(w, gconj(pi)); assert w is not None
        if ok:
            assert w in UNIT; data.append((z, a, w))
    for q in [q for q in range(3, 61, 2) if is_prime(q) and N % q]:
        ci = 1 if q % 8 in (1, 7) else -1      # (2/q)
        for (z0, a0, e0), (z1, a1, e1) in itertools.combinations(data, 2):
            ratio = gmulm(z1, ginvm(z0, q), q)
            k = (UNIT[e1] - UNIT[e0]) % 4
            pred = ci if k % 2 else 1
            for p, x1, x0 in zip(primes, a1, a0):
                if x1 != x0 and legendre(p, q) == -1: pred = -pred
            assert cls(ratio, q) == pred, (N, q, z0, z1)
            total += 1
    print("N=%d: %d primitive points; parity law verified for all pairs and all odd q<=60 coprime to N" % (N, len(data)))
print("R5: %d (pair, q) instances, all consistent with class(z_x/z_0) = class(eps_x/eps_0) prod chi_q(p_j)^(a_xj-a_0j)" % total)
