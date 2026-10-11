"""Referee-2 check T4.

A. Lemma 1.2 (validity of (R)) fails when v(c) = 0: exact counterexample on x^2+y^2 = 65.
B. "Twisted" point collisions (sum c = 2k != 0), a valid point-level residue input NOT included in Definition 1.2:
   for actual points,  prod_x z_x^(c_x) = E N^k A_v / conj(A_v),  v = sum c_x a_x - k e,  E = prod eps_x^(c_x).
   If prod rho_x^(c_x) = u N^k (mod q^a), q odd, q not dividing N, then q^a | E A_v - u conj(A_v).
   We check, exactly on actual points: the identity Z_+ conj(A) N^max(0,-k) = E A Z_- N^max(0,k), and
   (1+i) D'_u(c) | E A - u conj(A) != 0 whenever v != 0, for random c with sum c != 0 even.
C. Lemma T (twisted weights): in the profile P(H, pi, b, a*), every c with sum c = 2k != 0 has
   ||S^T c||_1 >= max(2 n, b (M - n)),  n = ||c||_1, and v = S^T c / 2 is outside L.  Exhaustive check on
   Paley-8 and Paley-12 profiles (b = 8) for entries in [-2,2], n <= 8 (M=8) / n <= 6 (M=12), plus c = lambda*1.
D. Completion with twisted residue constraints: for the parameters of Theorem 2 at P = P_* (q = 7, 11, 19, 43, 103,
   1003-ish) and for the five instances, the bad Haar proportion
       1.3201 (Pi - 1) + (8/pi)(1.0368/sqrt2) D_max sum_{n>=2 even} (2M)^n exp(-(tau/4) max(2n, b(M-n)))
   is < 1, where D_max = prod_{odd q <= Q0} q^(a_max) bounds every twisted collision modulus.
"""
import itertools, math, random
from decimal import Decimal, getcontext
from gz import *
getcontext().prec = 50
rng = random.Random(3)

def gcd_points_ok(pts, pis):
    return True

# ---------------- A ----------------
N = 65
pts = [(1, 8), (8, 1), (4, 7), (7, 4)]
assert all(norm(z) == N for z in pts)
lhs = mul(pts[0], pts[1]); rhs = mul(pts[2], pts[3])
print("A. (1+8i)(8+i) = %s, (4+7i)(7+4i) = %s" % (lhs, rhs))
pis = [gaussian_prime_over(5), gaussian_prime_over(13)]
def expo(z, pis):
    out = []
    for pi in pis:
        k = 0; w = z
        while divides(pi, w):
            w = div(w, pi); k += 1
        out.append(k)
    return out
rows = [expo(z, pis) for z in pts]
c = [1, 1, -1, -1]
v = [sum(ci * rw[j] for ci, rw in zip(c, rows)) for j in range(2)]
# collisions: prod rho^c = 1 mod every odd q^a <= 8 coprime to 65
Dv = 1
for q in (3, 7):
    a = 1
    while q ** (a + 1) <= 8:
        a += 1
    m = q ** a
    pr = (1, 0)
    for z, ci in zip(pts, c):
        pr = mulmod(pr, powmod(z, ci, m), m)
    if pr == (1, 0):
        Dv *= m
ang = [math.atan2(z[1], z[0]) for z in pts]
Cconst = (max(ang) - min(ang)) * N ** 0.5 / N ** 0.25
print("   exponent rows (pi_5, pi_13) = %s, c = %s, v(c) = %s, D_1(c) at Q0=8 = %d; (R) would need 0 >= %.2f; arc constant C = %.4f"
      % (rows, c, v, Dv, Dv / 2 ** 0.5, Cconst))
assert v == [0, 0] and Dv == 21

# ---------------- B ----------------
def circle_data(primes):
    N = 1
    for p in primes:
        N *= p
    pis = [gaussian_prime_over(p) for p in primes]
    data = []
    for z in lattice_points(N):
        a = expo(z, pis)
        base = (1, 0)
        for pi, k in zip(pis, a):
            base = mul(base, mul(power(pi, k), power(conj(pi), 1 - k)))
        eps = div(z, base)
        data.append((z, a, eps))
    return N, pis, data
units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
tot = 0; zero_v = 0
for primes in ((5, 13, 17, 29), (5, 13, 17, 29, 37), (13, 17, 29, 41, 53)):
    N, pis, data = circle_data(primes)
    Q0 = 60
    pps = []
    for q in primes_upto(Q0):
        if q == 2 or N % q == 0:
            continue
        a = 1
        while q ** (a + 1) <= Q0:
            a += 1
        pps.append((q, a))
    for _ in range(300):
        Msub = rng.randint(2, 6)
        sub = rng.sample(data, Msub)
        while True:
            cc = [rng.randint(-2, 2) for _ in range(Msub)]
            if any(cc) and sum(cc) % 2 == 0 and sum(cc) != 0:
                break
        k = sum(cc) // 2
        v = [sum(ci * d[1][j] for ci, d in zip(cc, sub)) - k for j in range(len(primes))]
        if not any(v):
            zero_v += 1
            continue
        A_ = (1, 0)
        for pi, vj in zip(pis, v):
            A_ = mul(A_, power(pi if vj > 0 else conj(pi), abs(vj)))
        E = (1, 0)
        for ci, d in zip(cc, sub):
            E = mul(E, power(d[2] if ci > 0 else conj(d[2]), abs(ci)))   # conj(unit) = unit^-1
        Zp = (1, 0); Zm = (1, 0)
        for ci, d in zip(cc, sub):
            if ci > 0:
                Zp = mul(Zp, power(d[0], ci))
            elif ci < 0:
                Zm = mul(Zm, power(d[0], -ci))
        Nk = (N ** abs(k), 0)
        left = mul(mul(Zp, conj(A_)), Nk if k < 0 else (1, 0))
        right = mul(mul(mul(E, A_), Zm), Nk if k > 0 else (1, 0))
        assert left == right
        for u in units:
            D = 1
            for (q, amax) in pps:
                best = 0
                for a in range(1, amax + 1):
                    m = q ** a
                    pr = (1, 0)
                    for ci, d in zip(cc, sub):
                        pr = mulmod(pr, powmod(d[0], ci, m), m)
                    targ = mulmod(u, (pow(N, k, m) if k >= 0 else pow(pow(N, -1, m), -k, m), 0), m)
                    if pr == targ:
                        best = a
                D *= q ** best
            X = (mul(E, A_)[0] - mul(u, conj(A_))[0], mul(E, A_)[1] - mul(u, conj(A_))[1])
            assert X != (0, 0)
            assert divides((1, 1), X) and X[0] % D == 0 and X[1] % D == 0
            assert norm(X) >= 2 * D * D
            tot += 1
print("B. twisted axiom (sum c != 0): identity and (1+i)D'|EA - u conj A != 0 verified on %d (character, unit) instances "
      "on 3 circles (skipped %d draws with v = 0)" % (tot, zero_v))

# ---------------- C ----------------
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
        H[I] = [-x for x in H[I]]
    return H
def profile(q, b):
    H = paley_norm(q); M = q + 1
    assign = [1] + list(range(1, M))
    cols = []
    for a in range(1, M):
        for _ in range(b):
            cols.append([H[x][a] for x in range(M)])
    for x in range(M):
        col = [H[y][assign[x]] for y in range(M)]; col[x] = -col[x]; cols.append(col)
    return cols
import numpy as np
for q, nmax in ((7, 8), (11, 6)):
    M = q + 1; b = 8
    S = np.array(profile(q, b), dtype=np.int64).T        # M x r
    worst = None; cnt = 0
    for k_ in range(1, nmax + 1):
        for supp in itertools.combinations(range(M), k_):
            for vv in itertools.product((-2, -1, 1, 2), repeat=k_):
                n = sum(map(abs, vv))
                if n > nmax or sum(vv) == 0 or sum(vv) % 2:
                    continue
                cvec = np.zeros(M, dtype=np.int64); cvec[list(supp)] = vv
                y1 = int(np.abs(cvec @ S).sum())
                bd = max(2 * n, b * max(0, M - n))
                sl = y1 - bd
                if worst is None or sl < worst[0]:
                    worst = (sl, n, list(vv))
                assert sl >= 0
                cnt += 1
    for lam in (1, 2, -1):
        cvec = lam * np.ones(M, dtype=np.int64)
        n = M * abs(lam)
        if (lam * M) % 2:
            continue
        assert int(np.abs(cvec @ S).sum()) >= max(2 * n, b * max(0, M - n))
    print("C. Paley-%d profile (b=8): ||S^T c||_1 >= max(2n, b(M-n)) on %d twisted c (entries in [-2,2], n<=%d) and c=lambda*1; min slack %d"
          % (M, cnt, nmax, worst[0]))

# ---------------- D ----------------
def first_primes_1mod4_above(x, cnt):
    out = []; t = x + 1
    while len(out) < cnt:
        if t % 4 == 1 and is_prime(t):
            out.append(t)
        t += 1
    return out
def Dmax(Q0):
    s = Decimal(0)
    for q in primes_upto(Q0):
        if q == 2:
            continue
        a = 1
        while q ** (a + 1) <= Q0:
            a += 1
        s += a * Decimal(q).ln()
    return s   # log D_max
def bad_total(M, b, tau, logDmax, Pi_minus_1):
    s = Decimal(0)
    n = 2
    while True:
        e = n * Decimal(2 * M).ln() - (tau / 4) * max(2 * n, b * max(0, M - n))
        term = e.exp()
        s += term
        if n > M and term < Decimal(10) ** -40:
            break
        n += 2
    tw = Decimal(8) / Decimal(math.pi) * Decimal("1.0368") / Decimal(2).sqrt() * logDmax.exp() * s
    return Decimal("1.3201") * Pi_minus_1 + tw, tw
print("D. completion bound with twisted residue constraints (must be < 1):")
for q in (7, 11, 19, 43, 103, 1019):
    M = q + 1; b = 8; r = b * (M - 1) + M
    L = first_primes_1mod4_above(2 * M, M)[-1]
    P = max(160 * M * M * L, 21 * r * r)
    tau = Decimal(P).ln()
    Pim1 = (Decimal(2 * r) / (Decimal(P).sqrt() - 1)).exp() - 1
    tot_, tw = bad_total(M, b, tau, Dmax(2 * M), Pim1)
    print("   Theorem 2 at P=P_*: M=%5d tau=%.3f  log D_max=%.1f  twisted part=%.3e  total=%.4f" % (M, tau, Dmax(2 * M), tw, tot_))
    assert tot_ < 1
for (q, P0exp) in ((11, 7), (19, 8), (23, 8), (43, 9), (11, 9)):
    M = q + 1; b = 8; r = b * (M - 1) + M
    ps = first_primes_1mod4_above(10 ** P0exp, r)
    tau = Decimal(ps[0]).ln()
    Pim1 = Decimal(1)
    for p in ps:
        Pim1 *= 1 + 2 / (Decimal(p).sqrt() - 1)
    Pim1 -= 1
    tot_, tw = bad_total(M, b, tau, Dmax(2 * M), Pim1)
    print("   instance M=%2d primes near 10^%d: twisted part=%.3e  total=%.4f" % (M, P0exp, tw, tot_))
    assert tot_ < 1
print("T4 DONE")
