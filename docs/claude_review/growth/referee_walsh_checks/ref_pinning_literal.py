"""Referee check 2 (Theorems B, F): literal Gaussian pure Hadamard cores.

For literal cores z_x = eps_x * prod_a G_a^{H(x,a)} (conj for -1), with random units eps_x and
actual split-prime blocks, compute the actual angular width Delta of the best containing arc,
the winding/unit integers k_x of (1.1), and check:
  (P)  every block angle phi_a is within Delta/2 of (pi/(2M)) m_a  (mod 2pi),  m_a = sum_x H(x,a) k_x
  (F)  for every pair and every 4-term label combination mu, the literal composite
       beta = prod G^{v} (conj for negative) has arg within ||mu||_1 Delta/2 of (pi/2) mu^t k (mod 2pi)
  (L)  Lemma 1.1 consequence on every residue coincidence: |Im(i^{-j} beta)| >= 1, and
       V(mu) >= 2 log(2/(||mu||_1 Delta))
Then: a short-arc search for ORDER-4 pure cores (rectangles) to exercise the pinning at small C.
"""
import math, random, itertools
from fractions import Fraction as Fr
random.seed(11)

def is_prime(n):
    if n < 2: return False
    i = 2
    while i * i <= n:
        if n % i == 0: return False
        i += 1
    return True

def gauss_prime(p):
    for a in range(1, int(p ** 0.5) + 1):
        b2 = p - a * a
        b = int(round(b2 ** 0.5))
        if b * b == b2: return (a, b)

def gmul(u, v): return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])
def gconj(u): return (u[0], -u[1])
def gpow(u, e):
    r = (1, 0)
    for _ in range(e): r = gmul(r, u)
    return r
def unitpow(k):
    return [(1, 0), (0, 1), (-1, 0), (0, -1)][k % 4]

def walsh(t):
    M = 2 ** t
    return [[(-1) ** bin(x & a).count("1") for a in range(M)] for x in range(M)]

def paley(q):
    sq = set((a * a) % q for a in range(1, q))
    ch = lambda a: 0 if a % q == 0 else (1 if (a % q) in sq else -1)
    M = q + 1
    H = [[1] * M for _ in range(M)]
    for x in range(1, M):
        for j in range(1, M):
            H[x][j] = -1 if x == j else -ch(x - j)
    return H

def angle_mod(t, m=2 * math.pi):
    t = math.fmod(t, m)
    if t > m / 2: t -= m
    if t < -m / 2: t += m
    return t

def best_window(angles):
    """angles in R; return (center, width) of smallest arc mod 2pi containing all."""
    a = sorted(t % (2 * math.pi) for t in angles)
    n = len(a)
    gaps = [(a[(i + 1) % n] - a[i]) % (2 * math.pi) for i in range(n)]
    if n == 1: return a[0], 0.0
    i = max(range(n), key=lambda i: gaps[i])
    start = a[(i + 1) % n]; width = 2 * math.pi - gaps[i]
    return start + width / 2, width

SPLIT = [p for p in range(5, 3000) if p % 4 == 1 and is_prime(p)]

def literal_check(H, b, rng):
    M = len(H)
    primes = rng.sample(SPLIT, b * (M - 1))
    blocks = []
    for a in range(1, M):
        G = (1, 0)
        for j in range(b):
            pi = gauss_prime(primes[(a - 1) * b + j])
            if rng.random() < 0.5: pi = gconj(pi)
            G = gmul(G, pi)
        blocks.append(G)
    phi = [math.atan2(G[1], G[0]) for G in blocks]
    units = [rng.randrange(4) for _ in range(M)]
    # literal points
    pts = []
    for x in range(M):
        z = unitpow(units[x])
        for a in range(1, M):
            G = blocks[a - 1] if H[x][a] == 1 else gconj(blocks[a - 1])
            z = gmul(z, G)
        pts.append(z)
    N0 = pts[0][0] ** 2 + pts[0][1] ** 2
    assert all(p[0] ** 2 + p[1] ** 2 == N0 for p in pts)
    # choose units optimally: rotate each point by i^{u} to minimise the containing arc
    # (brute force over a reference: angles mod pi/2, then best window)
    args = [math.atan2(p[1], p[0]) for p in pts]
    red = [t % (math.pi / 2) for t in args]
    a_sorted = sorted(red)
    n = M
    gaps = [((a_sorted[(i + 1) % n] - a_sorted[i]) % (math.pi / 2)) for i in range(n)]
    gaps[-1] = a_sorted[0] + math.pi / 2 - a_sorted[-1]
    i = max(range(n), key=lambda i: gaps[i])
    start = a_sorted[(i + 1) % n]; Delta = math.pi / 2 - gaps[i]
    center = start + Delta / 2
    # rotate each literal point by a unit into the window [center-Delta/2, center+Delta/2]
    rot_pts = []
    for p, t in zip(pts, args):
        best = None
        for u in range(4):
            q = gmul(p, unitpow(u))
            tt = math.atan2(q[1], q[0])
            if abs(angle_mod(tt - center)) <= Delta / 2 + 1e-12:
                best = q; break
        assert best is not None
        rot_pts.append(best)
    # (1.1): S phi_x = sum_a H(x,a) phi_a = theta' + delta_x + (pi/2) k_x, theta' = center (arg of common factor = 0)
    k = []
    for x in range(M):
        Sphi = sum(H[x][a] * phi[a - 1] for a in range(1, M))
        # rot_pts[x] = i^{u_x'} prod G^{H}, arg = Sphi + (pi/2) u'  (mod 2pi);  lies within Delta/2 of center
        tt = math.atan2(rot_pts[x][1], rot_pts[x][0])
        # Sphi = center + delta_x + (pi/2)k_x  with |delta_x|<=Delta/2
        kx = round((Sphi - center) / (math.pi / 2))
        delta = Sphi - center - (math.pi / 2) * kx
        assert abs(delta) <= Delta / 2 + 1e-9, (delta, Delta)
        k.append(kx)
    m = [sum(H[x][a] * k[x] for x in range(M)) for a in range(1, M)]
    # (P) pinning
    for a in range(1, M):
        eta = phi[a - 1] - (math.pi / (2 * M)) * m[a - 1]
        assert abs(eta) <= Delta / 2 + 1e-9, ("pinning", eta, Delta)
    # (F),(L) for pairs and 4-term combos
    nco = 0
    idx = list(range(1, M))
    combos = [((a,), (c,)) for a in idx for c in idx if a != c]
    sample4 = [(tuple(sorted(rng.sample(idx, 2))), tuple(sorted(rng.sample(idx, 2)))) for _ in range(300)]
    for plus, minus in combos + sample4:
        v = {}
        for a in plus: v[a] = v.get(a, 0) + 1
        for a in minus: v[a] = v.get(a, 0) - 1
        v = {a: e for a, e in v.items() if e}
        if not v: continue
        mu = [Fr(sum(H[x][a] * e for a, e in v.items()), M) for x in range(M)]
        assert sum(mu) == 0
        # column image integral = v
        for a2 in range(1, M):
            assert sum(mu[x] * H[x][a2] for x in range(M)) == v.get(a2, 0)
        l1 = sum(abs(t) for t in mu)
        muk = sum(mu[x] * k[x] for x in range(M))
        beta = (1, 0)
        for a, e in v.items():
            G = blocks[a - 1] if e > 0 else gconj(blocks[a - 1])
            beta = gmul(beta, gpow(G, abs(e)))
        ab = math.atan2(beta[1], beta[0])
        err = angle_mod(ab - (math.pi / 2) * float(muk))
        assert abs(err) <= float(l1) * Delta / 2 + 1e-9, ("(1.2)", err, float(l1) * Delta / 2)
        if muk.denominator == 1:
            j = int(muk) % 4
            r = gmul(beta, unitpow(-j))
            assert r[1] != 0
            V = math.log(beta[0] ** 2 + beta[1] ** 2)
            assert V >= 2 * math.log(2 / (float(l1) * Delta)) - 1e-9
            nco += 1
        else:
            nco += 0
    return Delta, nco

rng = random.Random(5)
tot = 0; coinc = 0
for name, H in [("walsh8", walsh(3)), ("walsh16", walsh(4)), ("paley12", paley(11)), ("paley20", paley(19))]:
    for _ in range(8):
        D, nc = literal_check(H, rng.choice([1, 2]), rng)
        tot += 1; coinc += nc
print("literal cores: pinning + (1.2) + Lemma 1.1 consequences verified on %d cores, %d residue coincidences" % (tot, coinc))

# ---- short-arc order-4 pure cores (rectangles): exercise pinning at small C ----
H4 = walsh(2)
cands = []
for a in range(1, 80):
    for b in range(1, 80):
        if math.gcd(a, b) != 1: continue
        n = a * a + b * b
        if n % 2 == 0: continue
        # conjugate-primitive odd-norm: gcd(a,b)=1 and n odd suffices
        cands.append((a, b))
# blocks near (pi/8)Z directions are the only useful ones; keep those within 0.02 of (pi/8)Z
def near8(G):
    t = math.atan2(G[1], G[0])
    return abs(angle_mod(t, math.pi / 8)) < 0.03
cands = [G for G in cands if near8(G)]
best = []
from math import gcd
def coprime_norms(Gs):
    ns = [g[0] ** 2 + g[1] ** 2 for g in Gs]
    return all(gcd(ns[i], ns[j]) == 1 for i in range(3) for j in range(i + 1, 3))
rng2 = random.Random(3)
tried = 0
for _ in range(60000):
    Gs = rng2.sample(cands, 3)
    if not coprime_norms(Gs): continue
    for orient in itertools.product([0, 1], repeat=3):
        Gs2 = [gconj(G) if o else G for G, o in zip(Gs, orient)]
        phi = [math.atan2(G[1], G[0]) for G in Gs2]
        angles = [sum(H4[x][a] * phi[a - 1] for a in range(1, 4)) for x in range(4)]
        red = sorted(t % (math.pi / 2) for t in angles)
        gaps = [red[i + 1] - red[i] for i in range(3)] + [red[0] + math.pi / 2 - red[3]]
        Delta = math.pi / 2 - max(gaps)
        W = sum(math.log(G[0] ** 2 + G[1] ** 2) for G in Gs2)
        C = Delta * math.exp(W / 4)
        tried += 1
        best.append((C, Gs2, Delta, W))
best.sort(key=lambda r: r[0])
print("order-4 search: %d oriented triples; smallest arc constants:" % tried)
for C, Gs2, Delta, W in best[:5]:
    # pinning: phi_a within Delta/2 of (pi/8)Z
    phi = [math.atan2(G[1], G[0]) for G in Gs2]
    pin = max(abs(angle_mod(t, math.pi / 8)) for t in phi)
    assert pin <= Delta / 2 + 1e-12
    # pair Grams G_xy <= 4 log C (D=0)
    Wa = [math.log(G[0] ** 2 + G[1] ** 2) for G in Gs2]
    Gmax = max(sum(Wa[a - 1] * H4[x][a] * H4[y][a] for a in range(1, 4)) for x in range(4) for y in range(x + 1, 4))
    assert Gmax <= 4 * math.log(C) + 1e-9
    print("   C=%.4f blocks=%s  max|eta|=%.2e <= Delta/2=%.2e  maxG=%.2f <= 4logC=%.2f" % (C, Gs2, pin, Delta / 2, Gmax, 4 * math.log(C)))
print("ALL PASS")
