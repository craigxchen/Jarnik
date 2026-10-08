# aux_checks.py -- exact checks for auxiliary.md.  Run: python3 aux_checks.py
# Every assertion below is an exact identity/inequality in Z or Z[i] on actual lattice points of
# actual circles x^2+y^2=N (N = product of split primes, possibly with exponents).
import math, random, itertools, sys
from fractions import Fraction
from aux_gauss import *

random.seed(20261008)
LOG = []
def report(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.append(s)

# ---------------------------------------------------------------------------------------------
# Test circles and clusters
CIRCLES = [
    Circle([5, 13, 17, 29, 37, 41, 53, 61, 73, 89]),
    Circle([5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 109]),
    Circle([5, 13, 17, 29, 37, 41, 53], [2, 1, 3, 1, 2, 1, 1]),
    Circle(split_primes(100, 200)[:10]),
    Circle([5, 13, 17, 29, 37, 41, 53, 61, 73], [1, 2, 1, 1, 1, 1, 2, 1, 1]),
]

def clusters(circ, M, how_many=3):
    bc = best_clusters(circ, M)
    return [(C, [w[0] for w in win], [w[1] for w in win], [w[2] for w in win]) for C, win in bc[:how_many]]

# ---------------------------------------------------------------------------------------------
# A. Rational-curve factorization of interpolation determinants:
#    det[f_l(z_j)] * 2^{mD} * prod z_j^D = det[G_l(z_j)] = V(z) * Phi(z), Phi in Z[i] symmetric,
#    and Phi(z) ~ Phi(z_1,...,z_1) = det[G_l^{(j)}(z_1)/j!]  (confluent / Wronskian value).
def poly_mul(P, Q):
    R = {}
    for a, ca in P.items():
        for b, cb in Q.items():
            R[a + b] = gadd(R.get(a + b, (0, 0)), gmul(ca, cb))
    return {k: v for k, v in R.items() if v != (0, 0)}

def poly_pow(P, n):
    R = {0: (1, 0)}
    for _ in range(n): R = poly_mul(R, P)
    return R

def G_of_monomial(a, b, D, N):
    # 2^D Z^D X^a Y^b with X=(Z^2+N)/(2Z), Y=-i(Z^2-N)/(2Z)
    P = poly_mul(poly_pow({2: (1, 0), 0: (N, 0)}, a), poly_pow({2: (1, 0), 0: (-N, 0)}, b))
    P = poly_mul(P, {D - a - b: (2 ** (D - a - b), 0)})
    unit = gpow((0, -1), b)
    return {k: gmul(v, unit) for k, v in P.items()}

def poly_eval(P, z):
    s = (0, 0)
    for k, c in P.items(): s = gadd(s, gmul(c, gpow(z, k)))
    return s

def poly_deriv_over_fact(P, j):
    # coefficients of P^{(j)}/j!
    return {k - j: (c[0] * math.comb(k, j), c[1] * math.comb(k, j)) for k, c in P.items() if k >= j}

def f_eval(f, z):
    x, y = z
    s = (0, 0)
    for (a, b), c in f.items():
        s = gadd(s, gmul(c, (x ** a * y ** b, 0)))
    return s

def check_A():
    report("== A. determinant factorization on the circle ==")
    count = 0; ratios = []
    for circ in CIRCLES[:4]:
        for m in (3, 4, 5):
            for C, zs, avecs, units in clusters(circ, m, 2):
                D = (m + 1) // 2 + random.randint(0, 1)
                mons = [(a, b) for a in range(D + 1) for b in range(D + 1 - a)]
                fs = []
                for l in range(m):
                    f = {}
                    for mon in random.sample(mons, min(len(mons), 4)):
                        f[mon] = (random.randint(-3, 3), random.randint(-3, 3))
                    fs.append(f)
                Gs = []
                for f in fs:
                    G = {}
                    for (a, b), c in f.items():
                        for k, v in G_of_monomial(a, b, D, circ.N).items():
                            G[k] = gadd(G.get(k, (0, 0)), gmul(c, v))
                    Gs.append(G)
                detf = det_gauss([[f_eval(f, z) for z in zs] for f in fs])
                detG = det_gauss([[poly_eval(G, z) for z in zs] for G in Gs])
                lhs = gmul(detf, (2 ** (m * D), 0))
                for z in zs: lhs = gmul(lhs, gpow(z, D))
                assert lhs == detG, "factor z^D 2^D identity failed"
                V = (1, 0)
                for i in range(m):
                    for j in range(i + 1, m): V = gmul(V, gsub(zs[j], zs[i]))
                if detG == (0, 0): continue
                Phi = gdivexact(detG, V)
                assert Phi is not None, "Vandermonde does not divide"
                Phi0 = det_gauss([[poly_eval(poly_deriv_over_fact(G, j), zs[0]) for j in range(m)] for G in Gs])
                if Phi0 != (0, 0):
                    ratios.append(math.sqrt(gnorm(gsub(Phi, Phi0)) / gnorm(Phi0)))
                count += 1
    report("A: %d determinants: 2^{mD} prod z^D det f = det G and V | det G exactly; "
           "max |Phi - Phi(z1,..,z1)|/|Phi(z1,..,z1)| = %.3g" % (count, max(ratios)))

# ---------------------------------------------------------------------------------------------
# B. Anchor-relative ('gap principle') identities
def check_B():
    report("== B. anchor-relative heights, Thue-type gap = pair inequality ==")
    n = 0
    for circ in CIRCLES:
        for m in (3, 4, 5):
            for C, zs, avecs, units in clusters(circ, m, 2):
                z0 = zs[0]; a0 = avecs[0]
                P = []; Q = []
                for z, a in zip(zs, avecs):
                    if z == z0: P.append(None); Q.append(None); continue
                    g = ggcd(z, z0)
                    Pi = gdivexact(z, g)
                    # z/z0 = unit * Pi / conj(Pi)
                    lhs = gmul(z, gconj(Pi)); rhs = gmul(z0, Pi)
                    assert any(lhs == gmul(u, rhs) for u in UNITS)
                    assert is_conj_primitive(Pi)
                    d0i = circ.dist(a, a0)
                    assert abs(math.log(gnorm(Pi)) - d0i) < 1e-9
                    P.append(Pi); Q.append(gnorm(Pi))
                for i in range(1, m):
                    for j in range(i + 1, m):
                        T = gmul(gconj(P[i]), P[j])
                        G = math.gcd(T[0], T[1])
                        # G = rational content; primitive part norm = e^{d_ij}
                        prim = (T[0] // G, T[1] // G)
                        assert abs(math.log(gnorm(prim)) - circ.dist(avecs[i], avecs[j])) < 1e-9
                        assert G * G * gnorm(prim) == Q[i] * Q[j]
                        # Im(conj(Pi)Pj) = G * t with |t| = |Im prim| >= 1 (prim not real: conj-primitive non-unit)
                        assert T[1] % G == 0 and prim[1] != 0
                        n += 1
    report("B: %d anchor pairs: P_i conj-primitive, N(P_i)=e^{d_0i}, G_ij^2 e^{d_ij} = Q_i Q_j, "
           "Im(conj(P_i)P_j) = G_ij t_ij with t_ij != 0" % n)

# ---------------------------------------------------------------------------------------------
# C. Perfect-power (Baker) certificates via Smith normal form of the difference lattice
def smith(A):
    """Smith normal form: returns (U, S, V) with U*A*V = S diagonal, U,V unimodular (lists of lists)."""
    m, n = len(A), len(A[0])
    S = [row[:] for row in A]
    U = [[int(i == j) for j in range(m)] for i in range(m)]
    V = [[int(i == j) for j in range(n)] for i in range(n)]
    def swap_rows(M_, i, j): M_[i], M_[j] = M_[j], M_[i]
    def swap_cols(M_, i, j):
        for row in M_: row[i], row[j] = row[j], row[i]
    t = 0
    while t < min(m, n):
        # find nonzero pivot with minimal abs value in submatrix
        piv = None
        for i in range(t, m):
            for j in range(t, n):
                if S[i][j] != 0 and (piv is None or abs(S[i][j]) < abs(S[piv[0]][piv[1]])):
                    piv = (i, j)
        if piv is None: break
        i, j = piv
        swap_rows(S, t, i); swap_rows(U, t, i); swap_cols(S, t, j); swap_cols(V, t, j)
        done = False
        while not done:
            done = True
            for i in range(t + 1, m):
                q = S[i][t] // S[t][t]
                if q:
                    S[i] = [x - q * y for x, y in zip(S[i], S[t])]; U[i] = [x - q * y for x, y in zip(U[i], U[t])]
                if S[i][t] != 0:
                    swap_rows(S, t, i); swap_rows(U, t, i); done = False; break
            if not done: continue
            for j in range(t + 1, n):
                q = S[t][j] // S[t][t]
                if q:
                    for row in S: row[j] -= q * row[t]
                    for row in V: row[j] -= q * row[t]
                if S[t][j] != 0:
                    swap_cols(S, t, j); swap_cols(V, t, j); done = False; break
            if not done: continue
            # divisibility condition
            for i in range(t + 1, m):
                for j in range(t + 1, n):
                    if S[i][j] % S[t][t] != 0:
                        S[t] = [x + y for x, y in zip(S[t], S[i])]; U[t] = [x + y for x, y in zip(U[t], U[i])]
                        done = False; break
                if not done: break
        if S[t][t] < 0:
            S[t] = [-x for x in S[t]]; U[t] = [-x for x in U[t]]
        t += 1
    return U, S, V

def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]

def inv_unimodular(V):
    # exact inverse of a unimodular integer matrix via Fractions
    n = len(V)
    A = [[Fraction(x) for x in row] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(V)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]; A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    out = [[A[i][n + j] for j in range(n)] for i in range(n)]
    assert all(x.denominator == 1 for row in out for x in row)
    return [[int(x) for x in row] for row in out]

def certificate(B, k):
    """minimal D>0 and integer y with y^T B = D e_k^T (B: rows = difference vectors), or None."""
    U, S, V = smith(B)
    m, n = len(B), len(B[0])
    c = V[k]  # e_k = (e_k V) V^{-1} = c V^{-1};  Lambda = { w V^{-1} : w in d_1 Z x ... x d_r Z x 0 }
    r = sum(1 for i in range(min(m, n)) if S[i][i] != 0)
    if any(c[i] != 0 for i in range(r, n)): return None
    D = 1
    for i in range(r):
        d = S[i][i]; D = D * (d // math.gcd(d, c[i])) // math.gcd(D, d // math.gcd(d, c[i]))
    x = [D * c[i] // S[i][i] for i in range(r)] + [0] * (m - r)
    # y^T B = x^T S V^{-1}... we have U B V = S  =>  B = U^{-1} S V^{-1};  y^T B = (y^T U^{-1}) S V^{-1}
    # want y^T U^{-1} = x^T  => y^T = x^T U
    y = [sum(x[i] * U[i][j] for i in range(m)) for j in range(m)]
    yB = [sum(y[i] * B[i][j] for i in range(m)) for j in range(n)]
    assert yB == [D * int(j == k) for j in range(n)], (yB, D, k)
    return D, y

def char_product(zs, yfull):
    """Pi_y = prod_{y>0} z^y prod_{y<0} conj(z)^{|y|}"""
    P = (1, 0)
    for z, yi in zip(zs, yfull):
        if yi > 0: P = gmul(P, gpow(z, yi))
        elif yi < 0: P = gmul(P, gpow(gconj(z), -yi))
    return P

def check_C():
    report("== C. perfect-power certificates (Smith normal form of the difference lattice) ==")
    stats = []
    for circ in CIRCLES:
        for m in (3, 4, 5, 6, 7):
            for C, zs, avecs, units in clusters(circ, m, 2):
                B = [[a - b for a, b in zip(avecs[i], avecs[0])] for i in range(1, m)]
                varying = [k for k in range(len(circ.primes)) if any(row[k] for row in B)]
                for k in varying:
                    res = certificate(B, k)
                    if res is None: stats.append(("notinspan", m)); continue
                    D, y = res
                    yfull = [-sum(y)] + y   # balanced: anchor coefficient
                    Pi = char_product(zs, yfull)
                    # remove rational content n; remainder must be unit * pi_k^D or unit*conj(pi_k)^D
                    nrat = math.gcd(Pi[0], Pi[1]); prim = (Pi[0] // nrat, Pi[1] // nrat)
                    pik = circ.pis[k]
                    ok = any(prim == gmul(u, gpow(pik, 2 * D)) or prim == gmul(u, gpow(gconj(pik), 2 * D)) for u in UNITS)
                    assert ok, "certificate character is not a D-th power"
                    stats.append((D, m, len(varying), sum(abs(t) for t in yfull)))
    full = [s for s in stats if s[0] != "notinspan"]
    report("C: %d certificates verified exactly (Pi_y = n * unit * pi_k^{2D}); "
           "%d prime coordinates not in Q-span of differences" % (len(full), len(stats) - len(full)))
    from collections import Counter
    report("C: distribution of D (multiplier) by cluster size m:",
           sorted(Counter((s[1], s[0]) for s in full).items())[:40])

# ---------------------------------------------------------------------------------------------
# D. Quadratic class statistics: chi_2(z/z_*) = chi_2(unit) * prod_k (p_k/q)^{a_k}
def legendre(a, q):
    a %= q
    if a == 0: return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1

def fq2_mul(x, y, q):  # F_q[i], q = 3 mod 4
    return ((x[0] * y[0] - x[1] * y[1]) % q, (x[0] * y[1] + x[1] * y[0]) % q)

def fq2_pow(x, n, q):
    r = (1, 0)
    while n:
        if n & 1: r = fq2_mul(r, x, q)
        x = fq2_mul(x, x, q); n >>= 1
    return r

def chi2_point(z, circ, q):
    """quadratic character of z/z_* in G_q, z_* = prod conj(pi_k)^{e_k}"""
    w = (1, 0)
    for pi, e in zip(circ.pis, circ.exps): w = gmul(w, gpow(pi, e))   # conj(z_*) = prod pi^e
    Ninv = pow(circ.N % q, q - 2, q)
    x = gmul(z, w); x = ((x[0] * Ninv) % q, (x[1] * Ninv) % q)        # = z / z_*  mod q
    if q % 4 == 3:
        r = fq2_pow(x, (q + 1) // 2, q)
        assert r in ((1, 0), (q - 1, 0))
        return 1 if r == (1, 0) else -1
    iota = next(t for t in range(2, q) if (t * t + 1) % q == 0)
    return legendre(x[0] + iota * x[1], q)

def check_D():
    report("== D. quadratic class statistics ==")
    n = 0
    for circ in CIRCLES:
        pts = circ.all_points()
        qs = [q for q in range(3, 120) if is_prime(q) and circ.N % q != 0]
        for q in qs:
            unit_chi = []
            for u in range(4):
                unit_chi.append(chi2_point(UNITS[u], Circle.__new__(Circle), q) if False else None)
            for z, avec, u in random.sample(pts, min(60, len(pts))):
                lhs = chi2_point(z, circ, q)
                prod = 1
                for p, a in zip(circ.primes, avec): prod *= legendre(p, q) ** a
                # unit character: chi_2(i^u)
                if q % 4 == 3:
                    ci = 1 if q % 8 == 7 else -1
                else:
                    ci = 1 if q % 8 == 1 else -1
                assert lhs == prod * ci ** u, (q, avec, u)
                n += 1
    report("D: %d (point, q) pairs: chi_2(z/z_*) = chi_2(i)^u prod_k (p_k/q)^{a_k} exactly" % n)

# ---------------------------------------------------------------------------------------------
# E. Ptolemy integer relation n_A a + n_B b = n_C c (two-valued model, squarefree N)
def check_E():
    report("== E. Ptolemy relation in the two-valued model ==")
    n = 0
    for circ in [c for c in CIRCLES if all(e == 1 for e in c.exps)]:
        for C, zs, avecs, units in clusters(circ, 4, 3) + clusters(circ, 5, 2):
            for quad in itertools.combinations(range(len(zs)), 4):
                z1, z2, z3, z4 = [zs[t] for t in quad]
                a1, a2, a3, a4 = [avecs[t] for t in quad]
                X = gmul(gsub(z1, z2), gsub(z3, z4)); Y = gmul(gsub(z1, z4), gsub(z2, z3)); Z = gmul(gsub(z1, z3), gsub(z2, z4))
                assert Z == gadd(X, Y)
                assert gmul(X, gconj(Z))[1] == 0 and gmul(Y, gconj(Z))[1] == 0   # collinear with 0
                nA = nB = nC = 1
                for k, p in enumerate(circ.primes):
                    s = (a1[k], a2[k], a3[k], a4[k])
                    if s[0] == s[1] and s[2] == s[3] and s[0] != s[2]: nA *= p
                    if s[0] == s[3] and s[1] == s[2] and s[0] != s[1]: nB *= p
                    if s[0] == s[2] and s[1] == s[3] and s[0] != s[1]: nC *= p
                # primitive direction omega of the common line
                g = math.gcd(math.gcd(X[0], X[1]), math.gcd(Y[0], Y[1]))
                omega = (Z[0] // math.gcd(Z[0], Z[1]), Z[1] // math.gcd(Z[0], Z[1]))
                def coef(T):
                    # T = t * omega, t integer
                    if omega[0] != 0:
                        assert T[0] % omega[0] == 0; t = T[0] // omega[0]
                    else:
                        t = T[1] // omega[1]
                    assert gmul((t, 0), omega) == T
                    return t
                x, yv, zc = coef(X), coef(Y), coef(Z)
                assert x % nA == 0 and yv % nB == 0 and zc % nC == 0
                assert zc == x + yv
                n += 1
    report("E: %d quadruples: X+Y=Z collinear with 0, n_A | x, n_B | y, n_C | z, so n_A a + n_B b = n_C c" % n)

# ---------------------------------------------------------------------------------------------
# F. Uniform weights are optimal for linear combinations of pair inequalities
def check_F():
    import numpy as np
    report("== F. optimal pair weights ==")
    worst = 0
    for trial in range(2000):
        M = random.randint(3, 12)
        rho = np.zeros((M, M))
        for i in range(M):
            for j in range(i + 1, M):
                v = random.random() ** random.choice([1, 3, 6]) if random.random() < 0.8 else 0.0
                rho[i, j] = rho[j, i] = v
        if rho.sum() == 0: continue
        mu = -np.linalg.eigvalsh(rho)[0]
        tot = rho.sum() / 2
        ratio = tot / (M * (M - 1) / 2 * mu)
        worst = max(worst, ratio)
        # exact max cut check for small M: maxcut - tot/2 <= M mu/4
        if M <= 10:
            best = 0
            for mask in range(1 << (M - 1)):
                x = np.array([1 if (mask >> t) & 1 else -1 for t in range(M)])
                best = max(best, -(x @ rho @ x) / 4)
            assert best <= M * mu / 4 + 1e-9
    assert worst <= 1 + 1e-9
    report("F: 2000 random weight matrices: |rho| <= C(M,2) mu(rho) (max ratio %.4f) and "
           "maxcut(rho)-|rho|/2 <= M mu(rho)/4" % worst)

def check_C2():
    report("== C2. certificate identity on arbitrary point sets with more points than varying primes ==")
    from collections import Counter
    cnt = Counter(); n = 0
    for primes, exps in [([5, 13, 17], [3, 2, 2]), ([5, 13, 17, 29], [2, 2, 1, 1]), ([13, 17, 29, 37], [1, 1, 1, 1]),
                         ([5, 29, 41], [4, 3, 2])]:
        circ = Circle(primes, exps)
        pts = circ.all_points()
        for trial in range(40):
            m = len(primes) + random.randint(1, 4)
            sel = random.sample(pts, m)
            zs = [t[0] for t in sel]; avecs = [t[1] for t in sel]
            B = [[a - b for a, b in zip(avecs[i], avecs[0])] for i in range(1, m)]
            for k in range(len(primes)):
                res = certificate(B, k)
                if res is None: continue
                D, y = res
                yfull = [-sum(y)] + y
                Pi = char_product(zs, yfull)
                nrat = math.gcd(Pi[0], Pi[1]); prim = (Pi[0] // nrat, Pi[1] // nrat)
                pik = circ.pis[k]
                assert any(prim == gmul(u, gpow(pik, 2 * D)) or prim == gmul(u, gpow(gconj(pik), 2 * D)) for u in UNITS)
                cnt[D] += 1; n += 1
    report("C2: %d certificates: primitive part of Pi_y is unit * pi_k^{+-2D} exactly; D distribution %s"
           % (n, sorted(cnt.items())[:15]))

def check_D2():
    report("== D2. Legendre-twisted collision bound on actual point sets ==")
    n = 0; tight = 0
    for circ in CIRCLES:
        pts = circ.all_points()
        for trial in range(30):
            M = random.randint(4, 40)
            sel = random.sample(pts, M)
            for q in [q for q in range(3, 60) if is_prime(q) and circ.N % q]:
                nu = q - 1 if q % 4 == 1 else q + 1
                classes = {}
                for z, a, u in sel:
                    key = (z[0] % q, z[1] % q)
                    classes[key] = classes.get(key, 0) + 1
                coll = sum(c * (c - 1) // 2 for c in classes.values())
                S = 0
                for z, a, u in sel:
                    v = legendre(2, q) ** u
                    for p, ak in zip(circ.primes, a): v *= legendre(p, q) ** ak
                    S += v
                assert len(classes) <= nu
                assert Fraction(coll) >= Fraction(M * M + S * S, 2 * nu) - Fraction(M, 2)
                if S * S > M: tight += 1
                n += 1
    report("D2: %d (set, q): collisions_q >= (M^2 + S_q^2)/(2 nu_q) - M/2 with S_q = sum_i (2/q)^{u_i} prod_k (p_k/q)^{a_ik}; "
           "%d cases with S_q^2 > M" % (n, tight))

def check_G():
    report("== G. Lemma D4(2) torsion divisibility and Lemma B height isolation ==")
    n = 0
    for p in split_primes(5, 400):
        pi = gaussian_prime_above(p)
        for l in [l for l in range(3, 80) if is_prime(l) and l != p]:
            # order of rho = pi/conj(pi) in (Z[i]/l)^*: smallest t with pi^t = conj(pi)^t mod l
            t = 1
            while True:
                X = gsub(gpow(pi, t), gpow(gconj(pi), t))
                if X[0] % l == 0 and X[1] % l == 0: break
                t += 1
            assert gnorm(X) % (l * l) == 0 and gnorm(X) <= 4 * p ** t
            n += 1
    m = 0
    for circ in [c for c in CIRCLES if all(e == 1 for e in c.exps)]:
        for M in (3, 4, 5, 6, 7):
            for C, zs, avecs, units in clusters(circ, M, 3):
                for a0 in avecs:
                    big = [a for a in avecs if a != a0 and circ.dist(a, a0) > 0.75 * circ.W + math.log(max(C, 1e-300))]
                    assert len(big) <= 1
                    m += 1
    report("G: %d (p,l): l^2 | Norm(pi^t - conj(pi)^t) <= 4 p^t at the order t of rho mod l; "
           "%d anchors: at most one point with d_0i > 3W/4 + log C" % (n, m))

def check_F2():
    report("== F2. Theorem D1 symmetrization inequalities (exact rationals, brute force) ==")
    n = 0
    for trial in range(300):
        M = random.randint(3, 8)
        rho = {}
        for i in range(M):
            for j in range(i + 1, M):
                rho[(i, j)] = Fraction(random.randint(0, 9)) if random.random() < 0.85 else Fraction(0)
        tot = sum(rho.values())
        if tot == 0: continue
        B = M * (M - 1) // 2
        # cut side: maxcut(rho) >= |rho| * floor(M/2)ceil(M/2) / binom(M,2)
        best = 0
        for mask in range(1 << M):
            cut = sum(v for (i, j), v in rho.items() if ((mask >> i) & 1) != ((mask >> j) & 1))
            best = max(best, cut)
        assert best * B >= tot * (M // 2) * ((M + 1) // 2)
        # sieve side: min over partitions into <= nu classes <= |rho| E(M,nu)/binom(M,2)
        for nu in (2, 3):
            q, r = divmod(M, nu)
            E = nu * q * (q - 1) // 2 + r * q
            mn = None
            for assign in itertools.product(range(nu), repeat=M):
                v = sum(w for (i, j), w in rho.items() if assign[i] == assign[j])
                if mn is None or v < mn: mn = v
            assert mn * B <= tot * E
        n += 1
    report("F2: %d random integer weightings: maxcut(rho) >= |rho| maxcut_u/binom(M,2) and "
           "min_partition coll_rho <= |rho| E(M,nu)/binom(M,2) (nu=2,3)" % n)

if __name__ == "__main__":
    check_A(); check_B(); check_C(); check_C2(); check_D(); check_D2(); check_E(); check_F(); check_F2(); check_G()
    report("ALL CHECKS PASSED")
