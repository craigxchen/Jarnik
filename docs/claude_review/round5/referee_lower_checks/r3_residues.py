"""Referee check r3: residue inputs of lower.md (exact Gaussian-integer arithmetic).

(A) Prop 1.2 (validity of (R)) on actual lattice points: for random subsets of the lattice points of
    x^2+y^2 = N (N odd, all primes 1 mod 4, with repeated primes), random zero-sum c with v(c) != 0:
      |U - u V|^2 * prod_j p_j^{|v_j|}  >=  2 m_u(c)^2 N^{n/2}
    where m_u is the collision modulus of prod z_x^{c_x} against the unit u at odd q <= X, q not | N
    (u = 1 is Def 1.1(e); u = i, -1, -i is the all-unit extension), and for pairs on a common level at
    p_j <= X the level residues r_x = z_x / (pi^a conj(pi)^(e-a)) are multiplied in.
    Also reports the v(c) = 0 cases (excluded by Def 1.1(e)) and whether (R) fails there.
(B) Lemma 4.1 (residue design) at M = 44, X = 88: norms, m_c <= m^gen, m^gen | Im(Omega^4)/4,
    Im(Omega^4) != 0, m^gen <= Y^(2n)/4; and the referee's 2-adic extension
    2^(a2-1) * m^gen_odd | Im(Omega^4)/4 where a2 = 2-adic collision level of prod rho^c against +-1.
"""
import random, math, itertools
from rl_common import isprime

random.seed(99)
FAIL = []

# ---------- Gaussian integers as (a, b) tuples of Python ints
def gmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])

def gconj(x):
    return (x[0], -x[1])

def gnorm(x):
    return x[0] * x[0] + x[1] * x[1]

def gpow(x, k):
    r = (1, 0)
    for _ in range(k):
        r = gmul(r, x)
    return r

def gdivexact(x, y):
    n = gnorm(y)
    t = gmul(x, gconj(y))
    assert t[0] % n == 0 and t[1] % n == 0
    return (t[0] // n, t[1] // n)

def gdivides(y, x):
    n = gnorm(y)
    t = gmul(x, gconj(y))
    return t[0] % n == 0 and t[1] % n == 0

def gmod(x, m):
    return (x[0] % m, x[1] % m)

def ginv_mod(x, m):
    n = gnorm(x) % m
    ni = pow(n, -1, m)
    c = gconj(x)
    return ((c[0] * ni) % m, (c[1] * ni) % m)

UNITS = [(1, 0), (0, 1), (-1, 0), (0, -1)]

def two_squares(p):
    for a in range(1, int(math.isqrt(p)) + 1):
        b2 = p - a * a
        b = math.isqrt(b2)
        if b * b == b2:
            return (a, b)
    raise ValueError

def odd_primes_upto(X):
    return [q for q in range(3, X + 1) if isprime(q)]

def Aq(X, q):
    a = 0
    while q ** (a + 1) <= X:
        a += 1
    return a

def coll_level(P, eta, q, A):
    """max a <= A with P == eta mod q^a"""
    a = 0
    while a < A:
        m = q ** (a + 1)
        if (P[0] - eta[0]) % m == 0 and (P[1] - eta[1]) % m == 0:
            a += 1
        else:
            break
    return a

# ================= (A) validity on actual lattice points
def lattice_points(N):
    pts = []
    a = 0
    while a * a <= N:
        b2 = N - a * a
        b = math.isqrt(b2)
        if b * b == b2:
            for s in ((a, b), (-a, b), (a, -b), (-a, -b), (b, a), (-b, a), (b, -a), (-b, -a)):
                pts.append(s)
        a += 1
    return sorted(set(pts))

def exps(z, primes):
    """exponent of pi_j in z for each j (pi_j fixed)"""
    out = []
    for (p, pi) in primes:
        k = 0
        w = z
        while gdivides(pi, w):
            w = gdivexact(w, pi); k += 1
        out.append(k)
    return out

checked = 0; v0_cases = 0; v0_fail = 0; unit_checked = 0; level_checked = 0
for (Nfac, X) in (([(5, 2), (13, 1), (17, 1), (29, 1)], 40), ([(5, 1), (13, 1), (17, 1), (29, 1), (37, 1)], 60),
                  ([(13, 2), (17, 1), (41, 1)], 60), ([(5, 3), (13, 2)], 50)):
    N = 1
    for p, e in Nfac:
        N *= p ** e
    primes = [(p, two_squares(p)) for p, e in Nfac]
    evec = [e for p, e in Nfac]
    pts = lattice_points(N)
    # keep only primitive-content points: gcd with N structure irrelevant; we use all points
    qs = [q for q in odd_primes_upto(X) if N % q != 0]
    for trial in range(1500):
        M = random.randint(2, 8)
        Z = random.sample(pts, M)
        A = [exps(z, primes) for z in Z]
        # c: random zero-sum, small entries
        c = [0] * M
        for _ in range(random.randint(1, 3)):
            x, y = random.sample(range(M), 2)
            c[x] += 1; c[y] -= 1
        if not any(c):
            continue
        n = sum(abs(t) for t in c)
        v = [sum(c[x] * A[x][j] for x in range(M)) for j in range(len(primes))]
        U = (1, 0); V = (1, 0)
        for x in range(M):
            if c[x] > 0:
                U = gmul(U, gpow(Z[x], c[x]))
            elif c[x] < 0:
                V = gmul(V, gpow(Z[x], -c[x]))
        assert gnorm(U) == gnorm(V) == N ** (n // 2)
        wexp = 1
        for j, (p, pi) in enumerate(primes):
            wexp *= p ** abs(v[j])
        for u in UNITS:
            # collision modulus m_u at odd q not dividing N (u = 1 is Def 1.1(e))
            m = 1
            for q in qs:
                Aq_ = Aq(X, q); mod = q ** Aq_
                P = (1, 0)
                for x in range(M):
                    zx = gmod(Z[x], mod)
                    if c[x] > 0:
                        P = gmod(gmul(P, gpow(zx, c[x])), mod)
                    elif c[x] < 0:
                        P = gmod(gmul(P, gpow(ginv_mod(zx, mod), -c[x])), mod)
                m *= q ** coll_level(P, u, q, Aq_)
            # level residues at p_j <= X for same-level pairs
            supp = [x for x in range(M) if c[x] != 0]
            if len(supp) == 2 and n == 2:
                x, y = supp
                if c[x] < 0:          # orient as U/V: x carries c = +1
                    x, y = y, x
                for j, (p, pi) in enumerate(primes):
                    if p <= X and A[x][j] == A[y][j]:
                        e = evec[j]; a = A[x][j]
                        Dj = gmul(gpow(pi, a), gpow(gconj(pi), e - a))
                        rx = gdivexact(Z[x], Dj); ry = gdivexact(Z[y], Dj)
                        Ap = Aq(X, p); mod = p ** Ap
                        # pair collision rx == u * ry  (for u=1 this is Def 1.1(e))
                        lev = coll_level(gmod(gmul(rx, ginv_mod(ry, mod)), mod), u, p, Ap)
                        m *= p ** lev
                        level_checked += 1
            diff = (U[0] - (u[0] * V[0] - u[1] * V[1]), U[1] - (u[0] * V[1] + u[1] * V[0]))
            lhs = gnorm(diff) * wexp
            rhs = 2 * m * m * N ** (n // 2)
            if any(v):
                checked += (u == (1, 0)); unit_checked += 1
                if lhs < rhs:
                    FAIL.append(("validity", N, Z, c, u, m))
            elif u == (1, 0):
                v0_cases += 1
                if lhs < rhs:
                    v0_fail += 1
print("(A) validity of (R): %d (point set, c) instances with v(c) != 0 (unit 1), %d with all units, "
      "%d level-residue factors used; v(c)=0 instances: %d, of which (R) fails: %d"
      % (checked, unit_checked, level_checked, v0_cases, v0_fail),
      "ok" if not any(f[0] == "validity" for f in FAIL) else "FAIL")

# ================= (B) Lemma 4.1 design
M = 44; X = 88
qs = odd_primes_upto(X)
ells = []
t = X + 1
while len(ells) < M:
    if isprime(t) and t % 4 == 1:
        ells.append(t)
    t += 1
Y = ells[-1]
varpi = [two_squares(l) for l in ells]
# N: product of 30 primes = 1 mod 4 above 10^6 (only its residues matter)
Nlist = []
t = 10 ** 6
while len(Nlist) < 30:
    if isprime(t) and t % 4 == 1:
        Nlist.append(t)
    t += 1
N = 1
for p in Nlist:
    N *= p
beta = {}
for q in qs:
    A = Aq(X, q); mod = q ** A
    found = None
    for a in range(mod):
        for b in range(mod):
            if (a * a + b * b - N) % mod == 0:
                found = (a, b); break
        if found:
            break
    assert found
    beta[q] = found
A2 = Aq(X, 2); mod2 = 2 ** A2
beta2 = None
for a in range(mod2):
    for b in range(mod2):
        if (a * a + b * b - N) % mod2 == 0:
            beta2 = (a, b); break
    if beta2:
        break
assert beta2

def rho(x, q, mod, bq):
    return gmod(gmul(bq, gmul(varpi[x], ginv_mod(gconj(varpi[x]), mod))), mod)

# norms
for q in qs:
    mod = q ** Aq(X, q)
    for x in range(M):
        if (gnorm(rho(x, q, mod, beta[q])) - N) % mod:
            FAIL.append(("norm", q, x))
ntests = 0; maxratio = 0.0; n2adic = 0
for trial in range(1500):
    kind = trial % 3
    c = [0] * M
    if kind == 0:
        x, y = random.sample(range(M), 2); c[x] = 1; c[y] = -1
    elif kind == 1:
        for _ in range(random.randint(2, 4)):
            x, y = random.sample(range(M), 2); c[x] += 1; c[y] -= 1
    else:
        for _ in range(random.randint(2, 6)):
            x, y = random.sample(range(M), 2); k = random.randint(1, 3); c[x] += k; c[y] -= k
    if not any(c):
        continue
    n = sum(abs(t) for t in c)
    Om = (1, 0)
    for x in range(M):
        if c[x] > 0:
            Om = gmul(Om, gpow(varpi[x], c[x]))
        elif c[x] < 0:
            Om = gmul(Om, gpow(gconj(varpi[x]), -c[x]))
    Om4 = gpow(Om, 4)
    I4 = Om4[1]
    if I4 == 0 or I4 % 4:
        FAIL.append(("Im", c)); continue
    target = abs(I4 // 4)
    m1 = 1; mg = 1
    for q in qs:
        A = Aq(X, q); mod = q ** A
        P = (1, 0)
        for x in range(M):
            r_ = rho(x, q, mod, beta[q])
            if c[x] > 0:
                P = gmod(gmul(P, gpow(r_, c[x])), mod)
            elif c[x] < 0:
                P = gmod(gmul(P, gpow(ginv_mod(r_, mod), -c[x])), mod)
        lv = [coll_level(P, u, q, A) for u in UNITS]
        if sum(1 for l in lv if l > 0) > 1:
            FAIL.append(("two units", q, c))
        m1 *= q ** lv[0]
        mg *= q ** max(lv)
    if m1 > mg or target % mg != 0 or 4 * mg > Y ** (2 * n):
        FAIL.append(("design", c, m1, mg))
    # 2-adic extension
    P = (1, 0)
    for x in range(M):
        r_ = gmod(gmul(beta2, gmul(varpi[x], ginv_mod(gconj(varpi[x]), mod2))), mod2)
        if c[x] > 0:
            P = gmod(gmul(P, gpow(r_, c[x])), mod2)
        elif c[x] < 0:
            P = gmod(gmul(P, gpow(ginv_mod(r_, mod2), -c[x])), mod2)
    a2 = max(coll_level(P, u, 2, A2) for u in UNITS)
    if a2 >= 1:
        n2adic += 1
        if target % (mg * 2 ** (a2 - 1)) != 0:
            FAIL.append(("2adic", c, a2))
    ntests += 1
    maxratio = max(maxratio, math.log(mg) / math.log(Y ** (2 * n) / 4))
print("(B) design M=44, X=88, Y=%d: %d characters; m_c <= m^gen | Im(Omega^4)/4 <= Y^(2n)/4; "
      "max log m^gen / log(Y^(2n)/4) = %.3f; 2-adic collisions seen: %d"
      % (Y, ntests, maxratio, n2adic), "ok" if not any(f[0] in ("design", "Im", "two units", "2adic", "norm") for f in FAIL) else "FAIL")
print("ALL OK" if not FAIL else "FAILURES: %s" % FAIL[:5])
