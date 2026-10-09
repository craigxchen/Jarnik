"""Referee check R4: construct.md Lemma 4.1 (residue design) and two scope questions.

(a) Lemma 4.1, independently: rho_x = s * omega_x / conj(omega_x) mod q^a.  For every zero-sum c with entries
    in {+-1,+-2}, ||c||_1 <= nmax, and every unit u, compute D_u(c) from the residue data and check
      * D_u(c) * (1+i) DIVIDES Lambda_c - u conj(Lambda_c) in Z[i]   (stronger than the norm inequality),
      * sqrt2 D_u <= |Lambda_c - u conj Lambda_c| <= 2 L^(n/2).
(b) 2-adic extension (not in construct's class, but in A.3's Q_M = lcm(1..2M)): with rho_x^(2^a) =
    s_2 omega_x/conj omega_x mod 2^a (Norm s_2 = N mod 2^a), collisions mod 2^a also divide
    Lambda - u conj Lambda, so |Lambda - u conj Lambda| >= 2^(a_2) D_odd and the same (5.1) suffices.
(c) Prime-level parity (construct §7.1): for an actual cluster, rho_x/rho_0 lies in the coset
    eps-class * prod_j chi_q(p_j)^(a_xj - a_0j) of K^2 in the norm-one group K mod q.  We test whether the
    construct's M=12 instance (primes from 10000121, ell_x = 29..137, k = 0) satisfies this, allowing any
    choice of units eps_x.  Exact integer arithmetic throughout."""
import itertools, math
import numpy as np

def is_prime(n):
    if n < 2: return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0: return n == p
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

def two_sq(p):
    for a in range(1, int(p**0.5) + 1):
        b2 = p - a*a; b = math.isqrt(b2)
        if b*b == b2: return (a, b)

def gmul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def gmulm(a, b, m): return ((a[0]*b[0] - a[1]*b[1]) % m, (a[0]*b[1] + a[1]*b[0]) % m)
def gconj(a): return (a[0], -a[1])
def gnorm(a): return a[0]*a[0] + a[1]*a[1]
def gpow(a, k):
    r = (1, 0)
    for _ in range(k): r = gmul(r, a)
    return r
def gpowm(a, k, m):
    r = (1, 0); a = (a[0] % m, a[1] % m)
    while k:
        if k & 1: r = gmulm(r, a, m)
        a = gmulm(a, a, m); k >>= 1
    return r
def ginvm(a, m):
    n = gnorm(a) % m; ni = pow(n, -1, m)
    c = gconj(a); return ((c[0]*ni) % m, (c[1]*ni) % m)
def gdivides(d, z):
    n = gnorm(d); w = gmul(z, gconj(d)); return w[0] % n == 0 and w[1] % n == 0

UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]

def setup(M, Q0, Nprimes):
    ells = []; t = Q0 + 1
    while len(ells) < M:
        if t % 4 == 1 and is_prime(t): ells.append(t)
        t += 1
    om = [two_sq(l) for l in ells]
    qs = [q for q in range(2, Q0 + 1) if is_prime(q)]
    maxpow = {}
    for q in qs:
        m = q
        while m * q <= Q0: m *= q
        maxpow[q] = m
    smod = {}
    for q, m in maxpow.items():
        Nm = 1
        for p in Nprimes: Nm = Nm * p % m
        found = None
        for x0 in range(m):
            for y0 in range(m):
                if (x0*x0 + y0*y0 - Nm) % m == 0 and (x0*x0 + y0*y0) % q:
                    found = (x0, y0); break
            if found: break
        smod[q] = found
    return ells, om, maxpow, smod

def rho(x, q, m, om, smod):
    return gmulm(gmulm(smod[q], om[x], m), ginvm(gconj(om[x]), m), m)

def check_lemma41(M, Q0, nmax):
    Nprimes = []; t = 10**6 + 1
    while len(Nprimes) < 30:
        if t % 4 == 1 and is_prime(t): Nprimes.append(t)
        t += 1
    ells, om, maxpow, smod = setup(M, Q0, Nprimes)
    L = max(ells)
    cnt = 0; maxD = 1; max2 = 0; worst_ratio = None
    vals = (-2, -1, 1, 2)
    for k in range(2, nmax + 1):
        for supp in itertools.combinations(range(M), k):
            for vv in itertools.product(vals, repeat=k):
                if sum(vv) or sum(abs(v) for v in vv) > nmax: continue
                n = sum(abs(v) for v in vv)
                lam = (1, 0)
                for x, v in zip(supp, vv):
                    lam = gmul(lam, gpow(om[x] if v > 0 else gconj(om[x]), abs(v)))
                for u in UNITS:
                    ucl = gmul(u, gconj(lam)); diff = (lam[0] - ucl[0], lam[1] - ucl[1])
                    assert diff != (0, 0)
                    Dodd = 1; a2 = 0
                    for q, mq in maxpow.items():
                        m = q; best = 1
                        while m <= mq:
                            prod = (1, 0)
                            for x, v in zip(supp, vv):
                                r_ = rho(x, q, m, om, smod)
                                prod = gmulm(prod, gpowm(r_ if v > 0 else ginvm(r_, m), abs(v), m), m)
                            if prod == (u[0] % m, u[1] % m): best = m
                            m *= q
                        if q == 2:
                            a2 = int(round(math.log2(best)))
                        else:
                            Dodd *= best
                    # divisibility in Z[i]
                    assert gdivides((Dodd, 0), diff)
                    assert gdivides((1, 1), diff)
                    assert gdivides((2**a2, 0), diff)
                    nd = gnorm(diff)
                    assert 2 * Dodd * Dodd <= nd <= 4 * L**n
                    assert (2**a2 * Dodd)**2 <= nd
                    maxD = max(maxD, Dodd); max2 = max(max2, a2)
                    cnt += 1
    return cnt, maxD, max2, L

for (M, nmax) in ((8, 6), (12, 6), (20, 4)):
    Q0 = 2 * M
    cnt, maxD, max2, L = check_lemma41(M, Q0, nmax)
    print("Lemma 4.1 (+2-adic): M=%d Q0=%d L=%d ||c||_1<=%d: %d (char,unit) instances; D_odd*(1+i) | diff and "
          "2^a2 | diff and 2 D^2 <= |diff|^2 <= 4 L^n all hold; max D_odd=%d, max a2=%d" % (M, Q0, L, nmax, cnt, maxD, max2))

# ---------------- (c) prime-level parity test on the construct's M=12 instance
def legendre(a, q):
    a %= q
    return 0 if a == 0 else (1 if pow(a, (q - 1)//2, q) == 1 else -1)

q_pal = 11; M = 12; b = 8
Hk = np.zeros((M, M), dtype=np.int64)
for i in range(M):
    for j in range(M):
        if i == j: Hk[i, j] = 1
        elif i == 0: Hk[i, j] = 1
        elif j == 0: Hk[i, j] = -1
        else: Hk[i, j] = legendre((j - 1) - (i - 1), q_pal)
H = Hk * Hk[:, [0]]
cols = []
for a in range(1, M):
    for _ in range(b): cols.append(H[:, a].copy())
for x in range(M):
    a = 1 if x == 0 else x
    col = H[:, a].copy(); col[x] = -col[x]; cols.append(col)
S = np.array(cols).T; r = S.shape[1]; A = (1 + S) // 2
ps = []; t = 10**7 + 1
while len(ps) < r:
    if t % 4 == 1 and is_prime(t): ps.append(t)
    t += 1
assert ps[0] == 10000121
ells = []; t = 25
while len(ells) < M:
    if t % 4 == 1 and is_prime(t): ells.append(t)
    t += 1
print("parity test, construct's M=12 instance: p_j = %d..%d (r=%d), ell_x = %s" % (ps[0], ps[-1], r, ells))
oddq = [q for q in range(3, 25) if is_prime(q)]
bad = 0; detail = []
# bit beta[q][x] = chi_q(l_x) chi_q(l_0) * prod_j chi_q(p_j)^(a_xj - a_0j)  (must equal the unit class of eps_x/eps_0)
unit_nontriv = {q: ((q % 8) in (3, 5)) for q in oddq}       # class of i in K/K^2 is -1 iff (2/q) = -1
feasible_rows = 0
for x in range(1, M):
    vec = {}
    for q in oddq:
        bit = legendre(ells[x], q) * legendre(ells[0], q)
        for j in range(r):
            e = int(A[x, j] - A[0, j])
            if e and legendre(ps[j], q) == -1: bit = -bit
        vec[q] = bit
    # unit options: eps_x/eps_0 in {1,-1} (class trivial) or {i,-i} (class -1 exactly at q with (2/q) = -1)
    ok_trivial = all(vec[q] == 1 for q in oddq)
    ok_i = all(vec[q] == (-1 if unit_nontriv[q] else 1) for q in oddq)
    feasible_rows += (ok_trivial or ok_i)
    detail.append((x, "".join('+' if vec[q] == 1 else '-' for q in oddq), ok_trivial or ok_i))
print("  odd q <= 24:", oddq, " (class of i nontrivial at q = 3,5 mod 8)")
for x, s, ok in detail:
    print("  row %2d: parity vector %s  consistent with some unit choice: %s" % (x, s, ok))
print("  rows consistent with prime-level residues (any units): %d of %d" % (feasible_rows, M - 1))
print("R4 DONE")
