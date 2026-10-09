"""Referee check (arith lens) of sharp.md Section 7 (negative claims).

S1  Prop. 7.1 (Siegel): for random integer labellings lambda (M = 6, 7; |lambda| <= B) there is
    c != 0, sum c = 0, sum c lambda = 0 with ||c||_inf <= (M max(1,B))^(2/(M-2)); exhaustive
    search over the last M-2 coordinates (the first two are then forced).
S2  Count exponent of {c in Z^M : ||c||_1 <= 4M}: (1/M) log # -> max_alpha [alpha log 2 + H(alpha)
    + 4 H(alpha/4)] (claimed 'about e^(5.3M)').
S3  The Poisson heuristic of 7.3: sum of 1/q over primes q in [e^y, e^(2y)] (claimed factor
    (log 2 / y) per collision; the Mertens value is log 2).
S4  Coherent design zeta_x = n_x + i, n_x^2 + 1 prime > X (free point residues as in lower.md
    Lemma 4.1, rho_x = zeta_x / conj zeta_x): exact collision moduli at all prime powers <= X.
    Shows that (a) the unit-1 pair moduli divide n_x - n_y (tiny), far below the 'bound'
    4 log Y used in 7.2; (b) the unit -1 pair moduli divide n_x n_y + 1 and are large; (c) support-4
    characters have large unit-1 moduli.  So 7.2 shows only that the a-priori bound is useless,
    not that every coherent design fails; and Prop. 7.4 needs its pair hypothesis for all units.
"""
import math, itertools, random
import numpy as np

rng = random.Random(11)


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


# ---------------- S1
def siegel_trial(M, B):
    lam = [rng.randint(-B, B) for _ in range(M)]
    while lam[0] == lam[1]:
        lam = [rng.randint(-B, B) for _ in range(M)]
    H = math.floor((M * max(1, max(abs(l) for l in lam))) ** (2 / (M - 2)) + 1e-12)
    rngH = np.arange(-H, H + 1)
    grids = np.meshgrid(*([rngH] * (M - 2)), indexing='ij')
    rest = np.stack([g.ravel() for g in grids], axis=1)            # c_2..c_{M-1}
    S0 = rest.sum(axis=1)
    S1 = rest @ np.array(lam[2:])
    # c0 + c1 = -S0 ; lam0 c0 + lam1 c1 = -S1  ->  c0 = (-S1 + lam1 S0)/(lam0 - lam1)
    den = lam[0] - lam[1]
    numr = -S1 + lam[1] * S0
    okint = numr % den == 0
    c0 = numr // den
    c1 = -S0 - c0
    good = okint & (np.abs(c0) <= H) & (np.abs(c1) <= H) & ((np.abs(rest).sum(axis=1) + np.abs(c0) + np.abs(c1)) > 0)
    return bool(good.any()), H


print("S1 Siegel (Prop. 7.1):")
for M, B in [(6, 1), (6, 10), (6, 100), (7, 10), (7, 60)]:
    res = [siegel_trial(M, B) for _ in range(20)]
    print(f"   M={M} B={B}: H={res[0][1]}, solution found in {sum(r[0] for r in res)}/20 random labellings")

# ---------------- S2
def H2(p):
    return 0.0 if p <= 0 or p >= 1 else -p * math.log(p) - (1 - p) * math.log(1 - p)


best = max((a * math.log(2) + H2(a) + 4 * H2(a / 4), a) for a in np.linspace(1e-6, 1 - 1e-6, 200001))
print(f"S2 (1/M) log #{{||c||_1 <= 4M}} -> {best[0]:.4f} (at alpha = {best[1]:.3f})")
for M in [20, 50, 100]:
    tot = sum(2 ** k * math.comb(M, k) * math.comb(4 * M, k) for k in range(M + 1))
    print(f"   exact M={M}: (1/M) log # = {math.log(tot)/M:.4f}")

# ---------------- S3
for y in [10, 12, 14]:
    lo, hi = math.exp(y), math.exp(2 * y) if y <= 10 else None
for y in [5, 7, 8]:
    lo, hi = int(math.exp(y)), int(math.exp(2 * y))
    s = np.ones(hi + 1, dtype=bool); s[:2] = False
    for i in range(2, int(hi ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    ps = np.nonzero(s)[0]
    ps = ps[(ps >= lo) & (ps <= hi)]
    print(f"S3 y={y}: sum_(e^y<=q<=e^(2y)) 1/q = {np.sum(1.0/ps):.4f}  (log 2 = {math.log(2):.4f}; log2/y = {math.log(2)/y:.4f})")

# ---------------- S4
X = 2000
K = 10 ** 5
cands = []
n = K
while len(cands) < 12:
    if is_prime(n * n + 1):
        cands.append(n)
    n += 1
print(f"S4 coherent design zeta_x = n_x + i, n_x^2+1 prime, n_x in [{cands[0]}, {cands[-1]}], X = {X}")
qs = [q for q in range(3, X + 1, 2) if is_prime(q)]


def smooth_part(v):
    """product of q^a over odd primes q <= X with q^a | v, q^a <= X (the collision modulus)"""
    v = abs(v)
    m = 1
    for q in qs:
        a = 0
        while v % q == 0 and q ** (a + 1) <= X:
            v //= q; a += 1
        m *= q ** a
    return m


def gm(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


# verify collision modulus = smooth part of D_eta by direct residue computation at a few q
def coll_mod(Omega, s):
    P, Qi = Omega
    D = [Qi, P, P - Qi, P + Qi][s]   # eta = 1, -1, i, -i   (lower.md Lemma 4.1 step 2)
    return smooth_part(D)


pair1, pairm1, pairi = [], [], []
for x, y in itertools.combinations(range(len(cands)), 2):
    Om = gm((cands[x], 1), (cands[y], -1))
    pair1.append(math.log(coll_mod(Om, 0)))
    pairm1.append(math.log(coll_mod(Om, 1)))
    pairi.append(max(math.log(coll_mod(Om, 2)), math.log(coll_mod(Om, 3))))
print(f"   pairs: log m (unit 1) max {max(pair1):.2f} [divides n_y - n_x, |.| <= {cands[-1]-cands[0]}];"
      f"  unit -1: mean {np.mean(pairm1):.2f} max {max(pairm1):.2f};  units +-i: mean {np.mean(pairi):.2f};"
      f"  a-priori bound 4 log Y = {4*math.log(cands[-1]**2+1):.1f};  log X = {math.log(X):.2f}")
# direct residue check of the unit-1 pair criterion at all prime powers <= X for one pair
x, y = 0, 5
good = True
for q in qs:
    a = 1
    while q ** a <= X:
        Q = q ** a
        zx, zy = (cands[x] % Q, 1), (cands[y] % Q, 1)
        # rho_x/rho_y = zx conj(zy) / (conj(zx) zy) == 1 mod Q  iff  Im(zx conj zy) = 0 mod Q
        num = gm(zx, (zy[0], -1)); den = gm((zx[0], -1), zy)
        lhs = ((num[0] - den[0]) % Q == 0) and ((num[1] - den[1]) % Q == 0)
        good &= lhs == ((cands[y] - cands[x]) % Q == 0)
        a += 1
print(f"   direct residue check of the unit-1 pair criterion q^a | n_y - n_x at all q^a <= X: {good}")
quad = []
for P4 in itertools.combinations(range(len(cands)), 4):
    a, b, c, d = P4
    Om = gm(gm((cands[a], 1), (cands[b], 1)), gm((cands[c], -1), (cands[d], -1)))  # c = e_a+e_b-e_c-e_d
    quad.append(math.log(coll_mod(Om, 0)))
print(f"   support-4 characters (unit 1): log m mean {np.mean(quad):.2f}, max {max(quad):.2f}, "
      f"fraction with log m > log X/2: {np.mean(np.array(quad) > math.log(X)/2):.2f}")
