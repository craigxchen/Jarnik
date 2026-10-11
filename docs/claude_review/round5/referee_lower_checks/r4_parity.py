"""Referee check r4: the prime-level parity law (lower.md Prop 8.1) and the scope of the barrier.

(A) Prop 8.1 on actual lattice points (exact): for odd q not dividing N and the norm-one group
    U = {u in (Z[i]/q)^* : u conj(u) = 1} (cyclic of order q - chi_4(q)), the class of z_x/z_0 in U/U^2
    equals class(eps_x/eps_0) + sum_j (a_xj - a_0j) nu_j(q), nu_j(q) = [p_j non-residue mod q],
    class(i) = [q = +-3 mod 8].  Classes computed as t^(|U|/2) in {+1,-1}.
(B) lower.md's design on a concrete canonical b = 8 instance (q = 11, M = 12, X = 2M = 24): point residues
    rho_x = beta varpi_x / conj(varpi_x), auxiliary primes the first 12 primes = 1 mod 4 above 24, profile
    primes p_j the first r = 88 primes = 1 mod 4 above 10^6.  For q = 7, 17, 23 (= +-1 mod 8, units are
    squares) the law forces class(rho_x/rho_0) = sum_j (a_xj - a_0j) nu_j(q).  Count rows satisfying it at
    all three q.  (Factored data are NOT part of Def 1.1; this only measures the scope point of Sec. 8.)
"""
import math, random
from rl_common import isprime, chi, paley_H
import numpy as np

random.seed(3)

def gmul(x, y, m):
    return ((x[0] * y[0] - x[1] * y[1]) % m, (x[0] * y[1] + x[1] * y[0]) % m)

def gpow(x, k, m):
    r = (1, 0); b = (x[0] % m, x[1] % m)
    while k:
        if k & 1:
            r = gmul(r, b, m)
        b = gmul(b, b, m); k >>= 1
    return r

def ginv(x, m):
    n = (x[0] * x[0] + x[1] * x[1]) % m
    ni = pow(n, -1, m)
    return ((x[0] * ni) % m, (-x[1] * ni) % m)

def cls(u, q):
    """class of a norm-one unit u mod q in U/U^2 as +-1"""
    order = q - (1 if q % 4 == 1 else -1)
    t = gpow(u, order // 2, q)
    assert t in ((1, 0), (q - 1, 0)), t
    return 1 if t == (1, 0) else -1

def two_squares(p):
    for a in range(1, math.isqrt(p) + 1):
        b = math.isqrt(p - a * a)
        if b * b == p - a * a:
            return (a, b)

def gdiv_exact(x, y):
    n = y[0] ** 2 + y[1] ** 2
    t = (x[0] * y[0] + x[1] * y[1], x[1] * y[0] - x[0] * y[1])
    if t[0] % n or t[1] % n:
        return None
    return (t[0] // n, t[1] // n)

# ---------------- (A)
FAIL = []
inst = 0
for Nfac in ([5, 13, 17, 29, 37], [5, 5, 13, 41], [13, 17, 17, 53]):
    N = 1
    for p in Nfac:
        N *= p
    ps = sorted(set(Nfac)); e = {p: Nfac.count(p) for p in ps}
    pis = {p: two_squares(p) for p in ps}
    pts = []
    for a in range(-math.isqrt(N), math.isqrt(N) + 1):
        b2 = N - a * a
        if b2 < 0:
            continue
        b = math.isqrt(b2)
        if b * b == b2:
            pts += [(a, b), (a, -b)] if b else [(a, 0)]
    def decompose(z):
        """z = eps * prod pi^a conj(pi)^(e-a); return (eps, a-vector)"""
        w = z; avec = []
        for p in ps:
            pi = pis[p]; k = 0
            while True:
                t = gdiv_exact(w, pi)
                if t is None:
                    break
                w = t; k += 1
            # remaining conj(pi) power
            cp = (pi[0], -pi[1]); kk = 0
            while True:
                t = gdiv_exact(w, cp)
                if t is None:
                    break
                w = t; kk += 1
            if k + kk != e[p]:
                return None
            avec.append(k)
        assert w in ((1, 0), (0, 1), (-1, 0), (0, -1))
        return w, avec
    dec = [(z, decompose(z)) for z in pts]
    dec = [(z, d) for z, d in dec if d is not None]   # primitive points (no rational factor p)
    z0, (eps0, a0) = dec[0]
    for q in [q for q in range(3, 60) if isprime(q) and N % q]:
        for z, (eps, a) in dec[1:]:
            u = gmul(z, ginv(z0, q), q)
            # make it norm one: z/z0 has norm 1 mod q already
            lhs = cls(u, q)
            k_eps = {(1, 0): 0, (0, 1): 1, (-1, 0): 2, (0, -1): 3}
            de = (k_eps[eps] - k_eps[eps0]) % 4
            unit_cls = cls(gpow((0, 1), de, q), q) if de else 1
            rhs = unit_cls
            for j, p in enumerate(ps):
                if (a[j] - a0[j]) % 2:
                    rhs *= chi(p, q)
            inst += 1
            if lhs != rhs:
                FAIL.append((N, q, z))
print("(A) parity law on actual primitive points: %d instances," % inst, "ok" if not FAIL else FAIL[:3])

# ---------------- (B)
q0 = 11; M = 12; b = 8; X = 2 * M
H = paley_H(q0)
alpha = list(range(q0)) + [0]
cols = []
for a in range(q0):
    rows = [x for x in range(M) if alpha[x] == a]
    for t in range(b):
        col = H[:, a + 1].copy()
        if t < len(rows):
            col[rows[t]] *= -1
        cols.append(col)
S = np.array(cols).T
A = (1 + S) // 2
r = S.shape[1]
P = []
t = 10 ** 6
while len(P) < r:
    if isprime(t) and t % 4 == 1:
        P.append(t)
    t += 1
ells = []
t = X + 1
while len(ells) < M:
    if isprime(t) and t % 4 == 1:
        ells.append(t)
    t += 1
qs = [q for q in range(3, X + 1) if isprime(q) and q % 8 in (1, 7)]
good_rows = 0
for x in range(1, M):
    ok = True
    for q in qs:
        lhs = chi(ells[x], q) * chi(ells[0], q)      # class of rho_x/rho_0 (beta cancels)
        rhs = 1
        for j in range(r):
            if (A[x, j] - A[0, j]) % 2:
                rhs *= chi(P[j], q)
        if lhs != rhs:
            ok = False
    good_rows += ok
print("(B) q = %s: rows x != 0 of the M=12 design consistent with the parity law at all these q: %d of %d"
      % (qs, good_rows, M - 1))
# how many of the 12 auxiliary-prime assignments could be fixed by permuting/choosing ells?
# requirement per row: chi_q(ell_x) = chi_q(ell_0) * prod_j chi_q(p_j)^(a_xj - a_0j) for q in qs
need = []
for x in range(M):
    v = []
    for q in qs:
        s = 1
        for j in range(r):
            if (A[x, j] - A[0, j]) % 2:
                s *= chi(P[j], q)
        v.append(s)
    need.append(tuple(v))
print("    required Legendre patterns of ell_x relative to ell_0 at q = %s:" % qs, sorted(set(need)))
print("ALL OK" if not FAIL else "FAIL")
