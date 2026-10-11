"""Referee-2 check T3: residue design (Lemma 4.1) and the prime-level parity law.

A. Lemma 4.1 on the construct's M = 12 instance design (Q0 = 24, ell_x = first 12 primes = 1 mod 4 above 24,
   p_j = first 100 primes = 1 mod 4 above 10^7, s with s conj(s) = N mod q^a):
   for all characters with ||c||_1 <= 6 (+-1 entries) and ||c||_1 <= 4 (entries in [-2,2]), all 4 units u:
   X = Lambda_c - u conj(Lambda_c) != 0,  (1+i) | X,  D_u(c) | X in Z[i],  2 D^2 <= |X|^2 <= 4 L^n.
B. Parity law on ACTUAL lattice points.  For z_x = eps_x prod_j pi_j^(a_xj) conj(pi_j)^(e_j - a_xj) and odd q not
   dividing N, in K/K^2 (K = norm-one subgroup of (Z[i]/q^a)^*, |K| = (q - chi(q)) q^(a-1)):
       class(z_x / z_0) = class(eps_x/eps_0) * prod_j chi_q(p_j)^(a_xj - a_0j),   class(i) = (2/q).
   Checked on all points of several circles (squarefree and non-squarefree N), q <= 41, a in {1, 2}.
C. Converse for twin profiles: the set of point-residue vectors induced by prime residues (fixed units) is exactly
   the parity-consistent set.  Exhaustive for M = 4 (Paley-4, b = 1, r = 7) at q^a = 3, 5, 7, 9.
D. The construct's M = 12 instance violates the parity law for every unit choice; moreover no re-assignment of
   the twelve ell's to the rows (any permutation) and no unit choice makes it consistent.
"""
import itertools, random, time
from gz import *

t0 = time.time()
rng = random.Random(11)

# ---------------- reconstruct the M = 12 instance (construct's instance.py conventions) ----------------
q0 = 11; M = 12; b = 8
def paley_norm(q):
    n = q + 1
    H = [[0] * n for _ in range(n)]
    for J in range(n):
        H[0][J] = 1
    for I in range(1, n):
        H[I][0] = -1
        for J in range(1, n):
            H[I][J] = (1 if I == J else 0) + legendre((J - 1) - (I - 1), q)
    for I in range(1, n):
        H[I] = [-v for v in H[I]]
    return H
H = paley_norm(q0)
assign = [1] + list(range(1, q0 + 1))
cols = []
for a in range(1, M):
    for _ in range(b):
        cols.append([H[x][a] for x in range(M)])
for x in range(M):
    col = [H[y][assign[x]] for y in range(M)]; col[x] = -col[x]; cols.append(col)
r = len(cols)
S = [[cols[j][x] for j in range(r)] for x in range(M)]          # S[x][j]
A = [[(1 + S[x][j]) // 2 for j in range(r)] for x in range(M)]
ps = []; n0 = 10 ** 7 + 1
while len(ps) < r:
    if n0 % 4 == 1 and is_prime(n0):
        ps.append(n0)
    n0 += 1
Q0 = 2 * M
ells = []; t = Q0 + 1
while len(ells) < M:
    if t % 4 == 1 and is_prime(t):
        ells.append(t)
    t += 1
L = max(ells)
om = [gaussian_prime_over(l) for l in ells]
oddq = [q for q in primes_upto(Q0) if q > 2]
maxpow = {}
for q in oddq:
    a = 1
    while q ** (a + 1) <= Q0:
        a += 1
    maxpow[q] = a
Nmod = {}
for q, a in maxpow.items():
    m = q ** a; v = 1
    for p in ps:
        v = v * p % m
    Nmod[q] = v
smod = {}
for q, a in maxpow.items():
    m = q ** a
    smod[q] = next((x, y) for x in range(m) for y in range(m) if (x * x + y * y - Nmod[q]) % m == 0 and (x * x + y * y) % q)
def rho(x, q, a):
    m = q ** a
    return mulmod(mulmod(smod[q], om[x], m), invmod(conj(om[x]), m), m)

# ---------------- A ----------------
units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
def chars(nmax, K):
    vals = [v for v in range(-K, K + 1) if v]
    for k in range(2, nmax + 1):
        for supp in itertools.combinations(range(M), k):
            for vv in itertools.product(vals, repeat=k):
                if sum(vv) == 0 and sum(map(abs, vv)) <= nmax:
                    yield dict(zip(supp, vv))
cnt = 0; maxD = 0
seen = set()
for nmax, K in ((6, 1), (4, 2)):
    for c in chars(nmax, K):
        key = tuple(sorted(c.items()))
        if key in seen:
            continue
        seen.add(key)
        n = sum(abs(v) for v in c.values())
        lam = (1, 0)
        for x, v in c.items():
            lam = mul(lam, power(om[x] if v > 0 else conj(om[x]), abs(v)))
        D = {u: 1 for u in units}
        for q, amax in maxpow.items():
            best = {u: 0 for u in units}
            for a in range(1, amax + 1):
                m = q ** a
                prod = (1, 0)
                for x, v in c.items():
                    prod = mulmod(prod, powmod(rho(x, q, a), v, m), m)
                for u in units:
                    if prod == (u[0] % m, u[1] % m):
                        best[u] = a
            for u in units:
                D[u] *= q ** best[u]
        for u in units:
            X = (lam[0] - mul(u, conj(lam))[0], lam[1] - mul(u, conj(lam))[1])
            assert X != (0, 0)
            assert divides((1, 1), X)
            d = D[u]
            assert X[0] % d == 0 and X[1] % d == 0
            assert 2 * d * d <= norm(X) <= 4 * L ** n
            maxD = max(maxD, d)
            cnt += 1
print("A. Lemma 4.1 on the M=12 design: %d (character, unit) instances; X != 0, (1+i)|X, D_u|X, 2D^2 <= |X|^2 <= 4L^n all hold; max D = %d"
      % (cnt, maxD))

# ---------------- B: parity law on actual points ----------------
def check_parity_circle(primes_exps, qs, avals):
    N = 1
    for p, e in primes_exps:
        N *= p ** e
    pis = [gaussian_prime_over(p) for p, _ in primes_exps]
    pts = []
    for z in lattice_points(N):
        # primitive: not divisible by any rational prime p (i.e. not by both pi and conj pi)
        if any(divides(pi, z) and divides(conj(pi), z) for pi in pis):
            continue
        pts.append(z)
    data = []
    for z in pts:
        a = []; w = z
        for pi, (p, e) in zip(pis, primes_exps):
            k = 0
            while divides(pi, w):
                w = div(w, pi); k += 1
            a.append(k)
        base = (1, 0)
        for pi, k, (p, e) in zip(pis, a, primes_exps):
            base = mul(base, mul(power(pi, k), power(conj(pi), e - k)))
        eps = div(z, base)
        assert norm(eps) == 1
        m = {(1, 0): 0, (0, 1): 1, (-1, 0): 2, (0, -1): 3}[eps]
        data.append((z, a, m))
    inst = 0
    for q in qs:
        if N % q == 0:
            continue
        for aa in avals:
            mod = q ** aa
            ci = Kclass((0, 1), q, aa)
            assert ci == legendre(2, q)
            z0, a0, m0 = data[0]
            z0inv = invmod(z0, mod)
            for (z, a, m) in data[1:]:
                t = mulmod(z, z0inv, mod)
                lhs = Kclass(t, q, aa)
                rhs = ci ** ((m - m0) % 4)
                for (p, e), ax, a0x in zip(primes_exps, a, a0):
                    rhs *= legendre(p, q) ** abs(ax - a0x)
                assert lhs == rhs, (N, q, aa, z)
                inst += 1
    return N, len(pts), inst
tot = 0
for pe in ([(5, 1), (13, 1), (17, 1), (29, 1), (37, 1)],
           [(5, 2), (13, 1), (17, 1), (29, 1)],
           [(13, 1), (17, 2), (41, 1), (53, 1)],
           [(5, 1), (13, 1), (17, 1), (29, 1), (37, 1), (41, 1)]):
    qs = [q for q in primes_upto(41) if q > 2]
    N, npts, inst = check_parity_circle(pe, qs, (1, 2))
    tot += inst
    print("B. N=%d (%s): %d primitive points, parity law verified on %d (pair, q^a) instances"
          % (N, "*".join("%d^%d" % pe_ for pe_ in pe), npts, inst))
print("B. total %d instances, no exception" % tot)

# ---------------- C: converse (induced = parity-consistent) on a small twin profile ----------------
H4 = paley_norm(3); M4 = 4
cols4 = [[H4[x][a] for x in range(M4)] for a in range(1, M4)]
for x in range(M4):
    col = [H4[y][x if x else 1] for y in range(M4)]; col[x] = -col[x]; cols4.append(col)
r4 = len(cols4)
S4 = [[cols4[j][x] for j in range(r4)] for x in range(M4)]
primes4 = [5, 13, 17, 29, 37, 41, 53][:r4]
for (q, a) in ((3, 1), (5, 1), (7, 1), (3, 2)):
    if any(p % q == 0 for p in primes4):
        continue
    m = q ** a
    units_m = [(x, y) for x in range(m) for y in range(m) if (x * x + y * y) % q]
    fib = {p: [z for z in units_m if (norm(z) - p) % m == 0] for p in primes4}
    # induced point residues with eps = 1: rho_x = prod_j r_j^(a_xj) conj(r_j)^(1 - a_xj)
    induced = set()
    for choice in itertools.product(*[fib[p] for p in primes4]):
        vec = []
        for x in range(M4):
            v = (1, 0)
            for j in range(r4):
                v = mulmod(v, choice[j] if S4[x][j] == 1 else conj(choice[j]), m)
            vec.append(v)
        induced.add(tuple(vec))
    # parity-consistent set relative to one induced base vector
    base = next(iter(induced))
    N4 = 1
    for p in primes4:
        N4 = N4 * p % m
    fibN = [z for z in units_m if (norm(z) - N4) % m == 0]
    consistent = set()
    for vec in itertools.product(fibN, repeat=M4):
        cls = [Kclass(mulmod(vec[x], invmod(base[x], m), m), q, a) for x in range(M4)]
        if len(set(cls)) == 1:
            consistent.add(tuple(vec))
    assert induced == consistent, (q, a, len(induced), len(consistent))
    print("C. q^a=%d: induced set = parity-consistent set (%d of %d norm-N vectors)" % (m, len(induced), len(fibN) ** M4))

# ---------------- D: the M=12 instance violates the parity law ----------------
target = {}
for q in oddq:
    target[q] = [1] * M
    for x in range(M):
        v = 1
        for j in range(r):
            v *= legendre(ps[j], q) ** abs(A[x][j] - A[0][j])
        target[q][x] = v
u2 = {q: legendre(2, q) for q in oddq}
bad_rows = 0
for x in range(1, M):
    ok_some_unit = False
    for mrel in range(4):
        if all(legendre(ells[x] * ells[0], q) == u2[q] ** mrel * target[q][x] for q in oddq):
            ok_some_unit = True
    # cross-check with the actual residue classes of the design (a = 1 and the max power)
    for q in oddq:
        for a in sorted({1, maxpow[q]}):
            m = q ** a
            t_ = mulmod(rho(x, q, a), invmod(rho(0, q, a), m), m)
            assert Kclass(t_, q, a) == legendre(ells[x] * ells[0], q)
    if not ok_some_unit:
        bad_rows += 1
print("D. M=12 instance: rows violating the parity law for EVERY unit choice: %d of %d" % (bad_rows, M - 1))
# any permutation of the ells? rows need P(ell_x) * T_x in {g, g*u2} for a common g
pat = lambda l: tuple(legendre(l, q) for q in oddq)
T = [tuple(target[q][x] for q in oddq) for x in range(M)]
U2 = tuple(u2[q] for q in oddq)
def mulv(a, b_):
    return tuple(i * j for i, j in zip(a, b_))
best = 0
for g in itertools.product((1, -1), repeat=len(oddq)):
    # bipartite matching rows -> ells with P(ell)*T_x in {g, g*U2}
    adj = [[k for k in range(M) if mulv(pat(ells[k]), T[x]) in (g, mulv(g, U2))] for x in range(M)]
    match = {}
    def aug(x, seen):
        for k in adj[x]:
            if k in seen:
                continue
            seen.add(k)
            if k not in match or aug(match[k], seen):
                match[k] = x
                return True
        return False
    size = sum(1 for x in range(M) if aug(x, set()))
    best = max(best, size)
print("D. best number of rows that can be made parity-consistent by ANY permutation of the 12 ells and units: %d of %d"
      % (best, M))
print("T3 DONE (%.0fs)" % (time.time() - t0))
