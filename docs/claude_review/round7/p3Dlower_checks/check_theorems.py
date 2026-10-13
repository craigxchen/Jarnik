"""Exact checks of Theorem 1.1, Theorem 2.1 and Theorem 4.3 of p3Dlower.md on small instances (level 1 and
prime powers q^2 for q <= 13). No proof depends on this; it guards against sign/convention slips."""
import random, itertools, math
from sympy import primerange, isprime, legendre_symbol
from tgroup import Tq, gauss_prime
random.seed(7)
ok = True
# ---------- Theorem 1.1 (1): q = +-3 mod 8, adding iota*u permutes units by <v,u> ----------
for q in [3, 5, 11, 13, 19, 29, 37, 43, 53, 59, 61, 67]:
    T = Tq(q); m = T.m; io = T.dlog_i()
    assert q % 8 in (3, 5) and io % 2 == 1
    r = 6
    for trial in range(200):
        l1 = [random.randrange(m) for _ in range(r)]
        nu = [random.randrange(2) for _ in range(r)]
        u = [(l1[j] - nu[j]) % 2 for j in range(r)]
        l = [(l1[j] + io * u[j]) % m for j in range(r)]
        assert all(l[j] % 2 == nu[j] for j in range(r))
        v = [random.randint(-3, 3) for _ in range(r)]
        A = sum(a * b for a, b in zip(v, l)) % m; B = sum(a * b for a, b in zip(v, l1)) % m
        sh = sum(a * b for a, b in zip(v, u))
        for s in range(4):
            if (A == (s * io) % m) != (B == ((s - sh) * io) % m): ok = False
print("Thm 1.1(1) units permuted:", ok)
# ---------- Theorem 2.1: twin decomposition identity ----------
M, r = 5, 14
J = list(range(2 * M, r)); f = list(range(M)); s_ = list(range(M, 2 * M)); eps = [random.choice((1, -1)) for _ in range(M)]
a = [[0] * r for _ in range(M)]
for y in range(M):
    for x in range(M):
        sig = random.randrange(2)
        a[y][s_[x]] = sig
        a[y][f[x]] = sig + (eps[x] if x == y else 0)
    for j in J: a[y][j] = random.randrange(2)
# ensure 0/1 entries: choose sig so that sig+eps in {0,1}
for y in range(M):
    for x in range(M):
        if x == y:
            a[y][s_[x]] = 0 if eps[x] == 1 else 1; a[y][f[x]] = a[y][s_[x]] + eps[x]
assert all(v in (0, 1) for row in a for v in row)
for y in range(M):
    for x in range(M):
        assert a[y][f[x]] - a[y][s_[x]] == (eps[x] if x == y else 0)
m = 1000003 * 4
t21 = True
for trial in range(300):
    ell = [random.randrange(m) for _ in range(r)]
    U = [(ell[f[x]] + ell[s_[x]]) % m for x in range(M)]; Lam = [(eps[x] * ell[f[x]]) % m for x in range(M)]
    for y in range(M):
        F = (sum(a[y][j] * ell[j] for j in J) + sum(a[y][s_[x]] * U[x] for x in range(M))) % m
        if (F + Lam[y]) % m != sum(a[y][j] * ell[j] for j in range(r)) % m: t21 = False
    v = [random.randint(-3, 3) for _ in range(r)]
    cp = [eps[x] * (v[f[x]] - v[s_[x]]) for x in range(M)]
    rhs = (sum(v[j] * ell[j] for j in J) + sum(v[s_[x]] * U[x] for x in range(M)) + sum(cp[x] * Lam[x] for x in range(M))) % m
    if rhs != sum(v[j] * ell[j] for j in range(r)) % m: t21 = False
print("Thm 2.1 (kappa = F + Lambda, identity (2.2)):", t21)
# ---------- Theorem 4.3: twisted coherent labels; collisions <=> q | D_s(W_{t,d}) ----------
def gmul(x, y): return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def conj(x): return (x[0], -x[1])
def Ds(W, s):
    X_, Y_ = W
    return [Y_, X_ - Y_, X_, X_ + Y_][s]
Xmax = 200
splits = [p for p in primerange(5, 5000) if p % 4 == 1]
ells = [2] + splits                                     # twist list: 2, 5, 13, ...
zp = [p for p in splits if p > Xmax][:5]                # zeta norms > X
zeta = [gauss_prime(p) for p in zp]
P = []
while len(P) < 5:
    c = random.randrange(10**6, 2 * 10**6)
    if c % 4 == 1 and isprime(c) and c not in P: P.append(c)
t43 = True; ncoll = 0
for q in primerange(3, Xmax + 1):
    T = Tq(q); m = T.m; io = T.dlog_i()
    tq = next(i for i, l in enumerate(ells) if l % q != 0 and legendre_symbol(l % q, q) == -1)
    rho = gauss_prime(ells[tq]); lr = T.dlog_ratio(*rho)
    coh = [T.dlog_ratio(*z) for z in zeta]
    nu = [1 if legendre_symbol(p % q, q) == -1 else 0 for p in P]
    delta = [(coh[j] % 2) ^ nu[j] for j in range(5)]
    lab = [(coh[j] + delta[j] * lr) % m for j in range(5)]
    assert all(lab[j] % 2 == nu[j] for j in range(5))
    for v in itertools.product(range(-2, 3), repeat=5):
        if not any(v): continue
        val = sum(x * y for x, y in zip(v, lab)) % m
        d = sum(x * y for x, y in zip(v, delta))
        W = (1, 0)
        for j in range(5):
            for _ in range(abs(v[j])): W = gmul(W, zeta[j] if v[j] > 0 else conj(zeta[j]))
        for _ in range(abs(d)): W = gmul(W, rho if d > 0 else conj(rho))
        for s in range(4):
            coll = (val == (s * io) % m)
            div = (Ds(W, s) % q == 0)
            if Ds(W, s) == 0: t43 = False
            if coll != div: t43 = False
            ncoll += coll
print("Thm 4.3 (collision at q with unit i^s  <=>  q | D_s(Z_v rho^d), D_s != 0):", t43, " collisions tested:", ncoll)
print("ALL OK" if (ok and t21 and t43) else "FAILURE")
