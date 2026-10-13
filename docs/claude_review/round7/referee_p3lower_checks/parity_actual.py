#!/usr/bin/env python3
"""Referee check of Prop. B.3 (forced Legendre parity) and of the Sec. 4 description of H for actual
lattice points, in exact Gaussian-integer arithmetic.  Independent code.
Points: z_x = i^(u_x) prod_j pi_j^(a_xj) conj(pi_j)^(e_j - a_xj), rows a_x distinct, e_j in {1,2}.
For odd prime powers Q = q^a (q not dividing N):
 (i)  t = z_x/z_y mod Q has t^(m_Q/2) = (-1)^(sum_j (a_xj-a_yj) nu_j(q) + lambda_q (u_x-u_y));
 (ii) unit-1 collisions z_x = z_y mod Q join equal parities, and their number is >= E(n0,m/2)+E(n1,m/2);
 (iii) collision mod Q <=> Q | reduced chord (z_x - z_y)/g_xy <=> zx conj(zy)/|g|^2 = |B|^2 mod Q;
 (iv) Q | Im(zx conj(zy))/|g|^2 alone characterises collisions with unit +-1 (not unit 1)."""
import itertools, random, math
from sympy import primerange, legendre_symbol, factorint

random.seed(2026)

def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0], -a[1])
def gpow(a, n):
    r = (1, 0)
    for _ in range(n): r = gmul(r, a)
    return r
def gmod(a, Q): return (a[0] % Q, a[1] % Q)
def gpowmod(a, n, Q):
    r = (1, 0); b = gmod(a, Q)
    while n:
        if n & 1: r = gmod(gmul(r, b), Q)
        b = gmod(gmul(b, b), Q); n >>= 1
    return r
def norm(a): return a[0]*a[0]+a[1]*a[1]
def ggcd(a, b):
    # Euclid in Z[i]
    while b != (0, 0):
        n = norm(b)
        num = gmul(a, gconj(b))
        qq = (round(num[0] / n), round(num[1] / n))
        # exact rounding with integers
        qq = ((2*num[0] + n) // (2*n), (2*num[1] + n) // (2*n))
        r = (a[0] - gmul(qq, b)[0], a[1] - gmul(qq, b)[1])
        a, b = b, r
    return a
def gdiv_exact(a, b):
    n = norm(b); num = gmul(a, gconj(b))
    assert num[0] % n == 0 and num[1] % n == 0
    return (num[0] // n, num[1] // n)

def gaussian_prime(p):
    for x in range(1, p):
        if (x*x + 1) % p == 0:
            g = ggcd((p, 0), (x, 1))
            assert norm(g) == p
            return g
    raise ValueError

def E(n, Q):
    b, r = divmod(n, Q)
    return Q * b * (b - 1) // 2 + r * b

def mQ(q, a): return (q - (1 if q % 4 == 1 else -1)) * q ** (a - 1)

fails = 0; npairs = 0; ncoll = 0; nmod = 0; im_pm = 0
for trial in range(6):
    ps = random.sample([5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113], random.choice([5, 6, 7]))
    es = [random.choice([1, 1, 2]) for _ in ps]
    pis = [gaussian_prime(p) for p in ps]
    N = 1
    for p, e in zip(ps, es): N *= p ** e
    M = random.choice([30, 45, 60])
    rows = set()
    allrows = list(itertools.product(*[range(e + 1) for e in es]))
    rows = random.sample(allrows, min(M, len(allrows)))
    M = len(rows)
    us = [random.randrange(4) for _ in rows]
    Z = []
    for a, u in zip(rows, us):
        z = gpow((0, 1), u)
        for pi, aj, e in zip(pis, a, es):
            z = gmul(z, gmul(gpow(pi, aj), gpow(gconj(pi), e - aj)))
        assert norm(z) == N
        Z.append(z)
    for q in primerange(3, 60):
        if N % q == 0: continue
        nu = [1 if legendre_symbol(p, q) == -1 else 0 for p in ps]
        lam = 1 if q % 8 in (3, 5) else 0
        sig = [(sum(a[j]*nu[j] for j in range(len(ps))) + lam*u) % 2 for a, u in zip(rows, us)]
        a_lev = 1
        while q ** a_lev <= 400:
            Q = q ** a_lev; m = mQ(q, a_lev); nmod += 1
            coll = 0
            for x in range(M):
                for y in range(x + 1, M):
                    npairs += 1
                    zx, zy = Z[x], Z[y]
                    # t = zx * conj(zy) * N^{-1} mod Q  (= zx/zy since zy conj(zy) = N)
                    Ninv = pow(N, -1, Q)
                    t = gmod(gmul(gmul(zx, gconj(zy)), (Ninv, 0)), Q)
                    if (norm(t) - 1) % Q: fails += 1
                    tp = gpowmod(t, m // 2, Q)
                    expect = (1, 0) if (sig[x] - sig[y]) % 2 == 0 else (Q - 1, 0)
                    if tp != expect: fails += 1
                    c = gmod((zx[0]-zy[0], zx[1]-zy[1]), Q) == (0, 0)
                    if c:
                        coll += 1; ncoll += 1
                        if sig[x] != sig[y]: fails += 1
                    # reduced chord criterion
                    g = ggcd(zx, zy)
                    ch = gdiv_exact((zx[0]-zy[0], zx[1]-zy[1]), g)
                    c2 = gmod(ch, Q) == (0, 0)
                    w = gdiv_exact(gmul(zx, gconj(zy)), (norm(g), 0))
                    Bn = N // norm(g)  # |B|^2 = |zx/g|^2
                    c3 = (w[1] % Q == 0) and ((w[0] - Bn) % Q == 0)
                    c4 = (w[1] % Q == 0)
                    cpm = c or (gmod((zx[0]+zy[0], zx[1]+zy[1]), Q) == (0, 0))
                    if not (c == c2 == c3): fails += 1
                    if c4 != cpm: fails += 1
                    if c4 and not c: im_pm += 1
            n1 = sum(sig); n0 = M - n1
            if coll < E(n0, m // 2) + E(n1, m // 2): fails += 1
            a_lev += 1
print(f"moduli checked: {nmod}, pair-modulus instances: {npairs}, unit-1 collisions: {ncoll}")
print(f"instances with Q | Im(...) but no unit-1 collision (unit -1 collisions): {im_pm}")
print("ALL ACTUAL-POINT PARITY CHECKS PASSED" if fails == 0 else f"FAILURES {fails}")
