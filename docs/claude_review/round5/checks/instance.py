"""Concrete fake circles of Theorem 2 (construct.md), skew-paired fully flipped Paley profile.

Usage: python3 instance.py q b C Q0factor P0exp
  q = 3 mod 4 prime (M = q+1), b copies, arc constant C, residue moduli Q0 = Q0factor*M, primes near 10^P0exp.

What is EXACT here (integers / Fractions / 50-digit Decimal logs):
  * primes p_j = 1 mod 4, distinct, all inside one window [P, P e^eta], eta <= 1/(4r);  p_j >= 21 r^2;
  * profile, rank M, twins (saturation);
  * residue data rho_x = s * omega_x / conj(omega_x) mod q^a for all odd q^a <= Q0, with Norm(s) = N mod q^a;
  * for every character c with ||c||_1 <= nmax: the collision moduli D_u(c) computed from the residue data,
    the bound D_u(c) <= 2 L^(n/2), and D_1(c) | Im(Lambda_c);
  * the explicit union bound of Theorem 2 for these parameters (failure probability < 1).
ILLUSTRATION (not proof): one random delta, checked on all characters with ||c||_1 <= nmax; random monomials
outside the character span on one random realisation.
"""
import sys, itertools, random
from fractions import Fraction
from decimal import Decimal, getcontext
import numpy as np
from lib import *
from gauss5 import *

getcontext().prec = 50
q = int(sys.argv[1]); b = int(sys.argv[2]); C = Fraction(sys.argv[3]); Q0f = int(sys.argv[4]); P0exp = int(sys.argv[5])
nmax = int(sys.argv[6]) if len(sys.argv) > 6 else 6
M = q + 1
H = normalize(paley1(q))
assign = np.array([1] + list(range(1, q + 1)))
S, lab, frow = flipped_profile(H, b, assign)
r = S.shape[1]
print("Paley q=%d M=%d b=%d r=%d C=%s" % (q, M, b, r, C))

# ---- primes in a narrow window
P0 = 10 ** P0exp
assert P0 >= 21 * r * r
ps = []; n0 = P0 + 1
while len(ps) < r:
    if n0 % 4 == 1 and is_prime(n0):
        ps.append(n0)
    n0 += 1
w = [Decimal(p).ln() for p in ps]
tau = w[0]; eta = w[-1] - w[0]
print("primes: %d..%d, log-window eta = %.3e (need <= 1/(4r) = %.3e), tau = %.4f" % (ps[0], ps[-1], eta, 1 / (4 * r), tau))
assert eta <= Decimal(1) / (4 * r)
W = sum(w)
Pi = 1.0
for p in ps:
    Pi *= 1 + 2 / (p ** 0.5 - 1)
print("W = %.3f,  W/(M log M) = %.3f,  completion: 1.33(Pi-1) = %.4f < 1" % (W, float(W) / (M * np.log(M)), 1.33 * (Pi - 1)))

# ---- profile checks: rank and twins
assert np.linalg.matrix_rank(S.astype(float)) == M
for x in range(M):
    j = b * (M - 1) + x                    # flipped column of row x
    a = lab[j]
    sib = (a - 1) * b                      # first unflipped copy of label a
    d = S[:, j] - S[:, sib]
    assert (np.flatnonzero(d) == [x]).all()
print("rank M and twins: OK")

# ---- residue data
Q0 = Q0f * M
pps = odd_prime_powers_upto(Q0)
ells = []; t = Q0 + 1
while len(ells) < M:
    if t % 4 == 1 and is_prime(t):
        ells.append(t)
    t += 1
L = max(ells)
omegas = [two_squares(l) for l in ells]
Nmods = {}
for (qq, a, m) in pps:
    Nm = 1
    for p in ps:
        Nm = Nm * p % m
    Nmods[m] = Nm
# s with Norm(s) = N mod m for every prime power m (then reduce compatibly: choose mod the maximal power of each q)
smod = {}
maxpow = {}
for (qq, a, m) in pps:
    if maxpow.get(qq, (0, 0))[0] < a:
        maxpow[qq] = (a, m)
for qq, (a, m) in maxpow.items():
    found = None
    for x0 in range(m):
        for y0 in range(m):
            if (x0 * x0 + y0 * y0 - Nmods[m]) % m == 0 and (x0 * x0 + y0 * y0) % qq != 0:
                found = (x0, y0); break
        if found:
            break
    smod[qq] = found
def rho(x, m, qq):
    s = smod[qq]
    om = omegas[x]
    v = gmulmod(gmulmod(s, om, m), ginv_mod(gconj(om), m), m)
    return v
# validity: norms
for (qq, a, m) in pps:
    for x in range(M):
        z = rho(x, m, qq)
        assert (gnorm(z) - Nmods[m]) % m == 0
print("residue data: %d odd prime powers <= Q0=%d, ell_x in [%d, %d], all rho_x in the conic fibre Norm = N: OK" % (len(pps), Q0, ells[0], L, ))

units = {(1, 0): "1", (0, 1): "i", (-1, 0): "-1", (0, -1): "-i"}
def Dvals(c):
    """collision moduli D_u(c) from the residue data, for the four units"""
    D = {u: 1 for u in units}
    for qq, (amax, mmax) in maxpow.items():
        best = {u: 0 for u in units}
        for a in range(1, amax + 1):
            m = qq ** a
            prod = (1, 0)
            for x in range(M):
                if c[x]:
                    prod = gmulmod(prod, gpowmod(rho(x, m, qq), int(c[x]), m), m)
            for u in units:
                if ((prod[0] - u[0]) % m, (prod[1] - u[1]) % m) == (0, 0):
                    best[u] = a
        for u in units:
            D[u] *= qq ** best[u]
    return D

def Lambda(c):
    lam = (1, 0)
    for x in range(M):
        if c[x] > 0:
            lam = gmul(lam, gpow(omegas[x], int(c[x])))
        elif c[x] < 0:
            lam = gmul(lam, gpow(gconj(omegas[x]), int(-c[x])))
    return lam

# ---- characters to test
def chars(nmax, K=2):
    vals = [v for v in range(-K, K + 1) if v]
    for k in range(2, nmax + 1):
        for supp in itertools.combinations(range(M), k):
            for vv in itertools.product(vals, repeat=k):
                if sum(vv) or sum(abs(v) for v in vv) > nmax:
                    continue
                c = np.zeros(M, dtype=np.int64); c[list(supp)] = vv
                yield c

rng = random.Random(2026)
u = [Fraction(rng.getrandbits(48), 1 << 48) for _ in range(M)]
lnL = Decimal(L).ln(); ln2 = Decimal(2).ln()
worst_axiom = None; worst_suff = None; maxD = 0; cnt = 0; worstD = None
for c in chars(nmax):
    n = int(np.abs(c).sum())
    y = S.T @ c
    slack = sum(wj * (abs(int(yj)) - 1) for wj, yj in zip(w, y)) / 2      # V_c - W/2, exact to 50 digits
    D = Dvals(c)
    lam = Lambda(c)
    assert lam[1] != 0
    assert abs(lam[1]) % D[(1, 0)] == 0                                   # D_1 | Im Lambda
    for uu, d in D.items():
        assert d <= 2 * L ** (n / 2) + 1e-9                                 # D_u <= 2 L^(n/2)
    cu = sum(Fraction(int(ci)) * ui for ci, ui in zip(c, u))
    assert cu != 0
    lhs = Decimal(abs(cu.numerator)).ln() - Decimal(cu.denominator).ln() - ln2      # log |X_c|/Delta
    # axiom (unit 1, the only non-trivial one):  |sin X_c| >= D_1 e^{-V_c/2}/2,  X_c = Delta*(c.u)/2, Delta = C e^{-W/4}
    rhs_ax = Decimal(D[(1, 0)]).ln() - ln2 - slack / 2 - Decimal(C.numerator).ln() + Decimal(C.denominator).ln()
    rhs_suff = (n * lnL) / 2 - slack / 2 - Decimal(C.numerator).ln() + Decimal(C.denominator).ln()
    m1 = lhs - rhs_ax - Decimal("0.001"); m2 = lhs - rhs_suff - Decimal("0.001")
    worst_axiom = m1 if worst_axiom is None or m1 < worst_axiom else worst_axiom
    worst_suff = m2 if worst_suff is None or m2 < worst_suff else worst_suff
    if D[(1, 0)] > maxD:
        maxD = D[(1, 0)]; worstD = (n, list(c))
    cnt += 1
print("characters with ||c||_1 <= %d, entries in [-2,2]: %d checked" % (nmax, cnt))
print("  max collision modulus D_1(c) = %d at n=%d (bound 2 L^(n/2) holds for all; D_1 | Im Lambda for all)" % (maxD, worstD[0]))
print("  ILLUSTRATION (one random delta): min log-margin, axiom (R): %.2f ; sufficient form with L^(n/2): %.2f" % (worst_axiom, worst_suff))
assert worst_axiom > 0 and worst_suff > 0

# ---- union bound of Theorem 2 (rigorous, with the proved G-bounds of Theorem 1)
def g_lower(n):
    if n == 2:
        return b - 4
    return min(2 * n + b - 8, (b - 3) * n / 2 + b)
tot = Decimal(0)
Cd = Decimal(C.numerator) / Decimal(C.denominator)
etaD = eta
for n in range(2, 4000, 2):
    gl = Decimal(g_lower(n))
    expo = (Decimal(n) * (Decimal(2 * M).ln() + lnL / 2) - (tau / 2) * gl / 2 + etaD * r / 4)
    # term: (2M)^n L^{n/2} * 4*1.0368*e^{-(V-W/2)/2}/C * (2 Delta n/pi + 2), Delta <= C e^{-W/4} <= 1e-100
    term = Decimal(4) * Decimal("1.0368") / Cd * expo.exp() * (Decimal(2) + Decimal(n) * Decimal("1e-100"))
    tot += term
# geometric tail n > 4000: for b = 8 and n >= 4, g(n) = 2n, so consecutive terms (n -> n+2) have ratio
# <= rho * (2 + (n+2)e-100)/(2 + n e-100) <= 1.0001 rho with rho = exp(2 log 2M + lnL - tau)
rho = (2 * Decimal(2 * M).ln() + lnL - tau).exp() * Decimal("1.0001")
assert b == 8 and rho < 1
tot += term * rho / (1 - rho)
print("union bound (Theorem 2, all characters, strengthened form):  P(failure) <= %.3e   (tail ratio %.2e)" % (tot, rho))
assert tot < 1

# ---- monomials outside the character span (illustration, generic realisation)
Sf = S.astype(float)
_, sv, Vt = np.linalg.svd(np.vstack([Sf, np.ones((1, r))]) if False else Sf)
# U = {u : S u in R 1}: null space of S plus one particular solution of S u = 1
null = Vt[M:].T
u1 = np.linalg.lstsq(Sf, np.ones(M), rcond=None)[0]
rng2 = np.random.default_rng(5)
phi = null @ rng2.uniform(-1e3, 1e3, null.shape[1]) + u1 * rng2.uniform(-1e3, 1e3)
logp = np.array([float(x) for x in w])
worst = 1e9
for _ in range(100000):
    k = rng2.integers(1, 6)
    v = np.zeros(r, dtype=np.int64)
    v[rng2.choice(r, size=k, replace=False)] = rng2.choice([-2, -1, 1, 2], size=k)
    ang = float(v @ phi) % (np.pi / 4)
    dist = min(ang, np.pi / 4 - ang)
    need = np.arcsin(np.exp(-float(np.abs(v) @ logp) / 2))
    worst = min(worst, dist / need)
print("ILLUSTRATION: 100000 random monomials (generic realisation): min dist/required = %.2f" % worst)
