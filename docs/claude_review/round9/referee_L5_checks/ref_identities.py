#!/usr/bin/env python3
"""Referee check (independent code) of L5.md section 2 identities, Lemma 3.1 and Theorem 3.2.

Twin profile built from scratch; labels parity-respecting with genuine Legendre parities of random
primes p_j = 1 mod 4. Group: Z/m_{q^a}, m_{q^a} = (q - chi(q)) q^(a-1), T = (m/4) Z/m.
Exact integer arithmetic only.
"""
import random, itertools, sys

random.seed(20261013)

def is_prime(n):
    if n < 2: return False
    for p in range(2, int(n**0.5) + 1):
        if n % p == 0: return False
    return True

def chi4(q): return 1 if q % 4 == 1 else -1
def mq(q, a=1): return (q - chi4(q)) * q**(a - 1)
def legendre(a, q):
    t = pow(a % q, (q - 1) // 2, q)
    return -1 if t == q - 1 else t

# --- random twin profile -------------------------------------------------
M, nJ = 6, 8
r = 2 * M + nJ
# columns: J = 0..nJ-1, f(x) = nJ + 2x, s(x) = nJ + 2x + 1
Jc = list(range(nJ)); fc = [nJ + 2 * x for x in range(M)]; sc = [nJ + 2 * x + 1 for x in range(M)]
eps = [random.choice([1, -1]) for _ in range(M)]
A = [[0] * r for _ in range(M)]
for x in range(M):
    for j in Jc: A[x][j] = random.randint(0, 1)
for z in range(M):
    col = [random.randint(0, 1) for _ in range(M)]
    # a_{.,f(z)} - a_{.,s(z)} = eps_z e_z  with 0/1 entries
    if eps[z] == 1: col[z] = 0
    else: col[z] = 1
    for y in range(M):
        A[y][sc[z]] = col[y]
        A[y][fc[z]] = col[y] + (eps[z] if y == z else 0)
for z in range(M):
    for y in range(M):
        assert A[y][fc[z]] - A[y][sc[z]] == (eps[z] if y == z else 0)
        assert A[y][fc[z]] in (0, 1)
n = r - M
sigma = [[A[y][sc[x]] for x in range(M)] for y in range(M)]
phi = [[A[y][j] for j in Jc] + sigma[y] for y in range(M)]   # in Z^n

primes = []
p = 5
while len(primes) < r:
    if is_prime(p) and p % 4 == 1: primes.append(p)
    p += 4 if random.random() < 0.5 else 8
    while p % 4 != 1: p += 1
qs = [q for q in range(3, 60) if is_prime(q) and q not in primes]

def free_part(ell, m):
    U = [(ell[fc[x]] + ell[sc[x]]) % m for x in range(M)]
    return [ell[j] % m for j in Jc] + U

def twin_coords(v):
    w = [v[j] for j in Jc] + [v[sc[x]] for x in range(M)]
    c = [eps[x] * (v[fc[x]] - v[sc[x]]) for x in range(M)]
    return w, c

def dot(a, b, m): return sum(x * y for x, y in zip(a, b)) % m

fails = 0; cnt21 = cnt23 = cnt22 = 0
for q in qs:
    for a in (1, 2):
        m = mq(q, a); iota = m // 4
        nu = [1 if legendre(p, q) == -1 else 0 for p in primes]
        for trial in range(20):
            ell = [(nu[j] + 2 * random.randrange(m // 2)) % m for j in range(r)]
            f = free_part(ell, m)
            Lam = [(eps[x] * ell[fc[x]]) % m for x in range(M)]
            F = [dot(phi[y], f, m) for y in range(M)]
            kap = [dot(A[y], ell, m) for y in range(M)]
            for y in range(M):
                if kap[y] != (F[y] + Lam[y]) % m: fails += 1
            # (2.2) for random v
            for _ in range(10):
                v = [random.randint(-2, 2) for _ in range(r)]
                w, c = twin_coords(v)
                if dot(v, ell, m) != (dot(w, f, m) + dot(c, Lam, m)) % m: fails += 1
            # Prop 2.1: (w, eps e_y) collides with unit s  <=>  kappa_y = <phi_y - eps w, f> + eps s iota
            for _ in range(10):
                y = random.randrange(M); e = random.choice([1, -1])
                w = [random.randint(-1, 1) if random.random() < 0.3 else 0 for _ in range(n)]
                for s in range(4):
                    lhs = ((dot(w, f, m) + e * Lam[y]) % m) == (s * iota) % m
                    rhs = kap[y] == (dot([pp - e * ww for pp, ww in zip(phi[y], w)], f, m) + e * s * iota) % m
                    if lhs != rhs: fails += 1
                    cnt21 += 1
            # Prop 2.3: force kappa_x = kappa_y by re-choosing Lambda_y (product structure, Prop 2.2)
            x, y = random.sample(range(M), 2)
            if (kap[x] - F[y] - nu[fc[y]]) % 2 == 0:   # same parity class needed
                Lam2 = Lam[:]; Lam2[y] = (kap[x] - F[y]) % m
                kap2 = [(F[u] + Lam2[u]) % m for u in range(M)]
                assert kap2[x] == kap2[y]
                for _ in range(10):
                    w = [random.randint(-1, 1) if random.random() < 0.3 else 0 for _ in range(n)]
                    e = random.choice([1, -1])
                    w2 = [ww + e * (px - py) for ww, px, py in zip(w, phi[x], phi[y])]
                    for s in range(4):
                        l = ((dot(w, f, m) + e * Lam2[y]) % m) == (s * iota) % m
                        rr = ((dot(w2, f, m) + e * Lam2[x]) % m) == (s * iota) % m
                        if l != rr: fails += 1
                        cnt23 += 1
            # Prop 2.2: ell recovered from (f, Lambda)
            ell_back = [0] * r
            for j_i, j in enumerate(Jc): ell_back[j] = f[j_i]
            for xx in range(M):
                ell_back[fc[xx]] = (eps[xx] * Lam[xx]) % m
                ell_back[sc[xx]] = (f[nJ + xx] - eps[xx] * Lam[xx]) % m
            if ell_back != ell: fails += 1
            cnt22 += 1

# Cor 2.4 norm inequality, all pairs
c24 = 0
for x in range(M):
    for y in range(M):
        if x == y: continue
        v = [A[x][j] - A[y][j] for j in range(r)]
        dphi = [a - b for a, b in zip(phi[x], phi[y])]
        if not (sum(map(abs, dphi)) >= (sum(map(abs, v)) - 4) / 2): fails += 1
        c24 += 1

print('Prop 2.1 instances', cnt21, '; Prop 2.3 instances', cnt23, '; Prop 2.2 round trips', cnt22,
      '; Cor 2.4 pairs', c24, '; failures', fails)

# Remark 2 of section 1: c' = +-e_y never the twin part of a character sum c_x a_x, sum c = 0
# (twin part of sum c_x a_x equals c): checked symbolically: twin_coords(a_x) = (phi_x, e_x)
for x in range(M):
    w, c = twin_coords(A[x])
    assert w == phi[x] and c == [1 if u == x else 0 for u in range(M)], 'twin coords of a_x'
print('twin coordinates of a_x are (phi_x, e_x): ok')
sys.exit(1 if fails else 0)
