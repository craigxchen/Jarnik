#!/usr/bin/env python3
"""Exact checks for Section 3 of hardness.md (single linear forms, products, rational directions).

(E) The 10-point cluster on N = 1176852625 (referee_families.md): for every pair, the pair norm
    n = N/|gcd|^2, its factorisation, the actual angle, the elementary (Liouville) lower bound
    sqrt2/sqrt(n), and the log10 of Matveev's lower bound for the same linear form.
(F) Product (cube) clusters: for random coprime block triples, the 8 product points satisfy
    C >= sqrt2 * N^(1/12) (Lemma 3.3: m = 3 blocks force N <= (C^2/2)^6).
(G) Rational directions: an arc within angle kappa/sqrt(R) of the direction of a Gaussian integer w
    holds at most 2|w|(kappa+C)^2 + 2 points (Lemma 4.2).
Integers are exact; floating point is used only for angles and logarithms in the printed tables.
"""
import math, random, itertools
random.seed(7)
out = []
def log(s):
    print(s); out.append(s)

def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0], -a[1])
def gnorm(a): return a[0]*a[0]+a[1]*a[1]
def gdivmod_round(a, b):
    # a / b rounded to nearest Gaussian integer
    nb = gnorm(b); num = gmul(a, gconj(b))
    q = ((2*num[0] + nb)//(2*nb), (2*num[1] + nb)//(2*nb))
    r = (a[0]-gmul(q,b)[0], a[1]-gmul(q,b)[1])
    return q, r
def ggcd(a, b):
    while b != (0, 0):
        _, r = gdivmod_round(a, b); a, b = b, r
    return a

def factor(n):
    f = {}; d = 2
    while d*d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0)+1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0)+1
    return f

def circle_points(N):
    pts = []
    x = 0
    r = math.isqrt(N)
    for x in range(-r, r+1):
        y2 = N - x*x
        y = math.isqrt(y2)
        if y*y == y2:
            pts.append((x, y))
            if y: pts.append((x, -y))
    return pts

def best_cluster(N, k):
    """smallest arc (in units of sqrt R) containing k lattice points of x^2+y^2=N."""
    pts = circle_points(N)
    R = math.sqrt(N)
    ang = sorted((math.atan2(p[1], p[0]) % (2*math.pi), p) for p in pts)
    m = len(ang)
    best = (1e300, None)
    for i in range(m):
        j = (i + k - 1)
        a0 = ang[i][0]; a1 = ang[j % m][0] + (2*math.pi if j >= m else 0)
        arc = R*(a1-a0)
        if arc < best[0]:
            best = (arc, [ang[(i+t) % m][1] for t in range(k)])
    return best[0]/math.sqrt(R), best[1]

split_primes = [p for p in range(5, 400) if p % 4 == 1 and all(p % q for q in range(2, int(p**0.5)+1))]
def gauss_prime(p):
    for a in range(1, math.isqrt(p)+1):
        b2 = p - a*a; b = math.isqrt(b2)
        if b*b == b2: return (a, b)

# (E)
N = 1176852625
C10, P = best_cluster(N, 10)
R = math.sqrt(N)
log(f"(E) N={N}={factor(N)}: 10-point cluster with C={C10:.6f}")
def matveev_log10(nlogs, A_prod, B, D=2):
    kappa = 2
    Cn = min((1/kappa)*(math.e*nlogs/2)**kappa*30**(nlogs+3)*nlogs**3.5, 2**(6*nlogs+20))
    lnlow = -Cn*D*D*(1+math.log(D))*(1+math.log(max(B, math.e)))*A_prod
    return lnlow/math.log(10)
rows = []
sqfree_pairs = 0
for z, w in itertools.combinations(P, 2):
    g = ggcd(z, w)
    npair = N // gnorm(g)
    assert N % gnorm(g) == 0
    f = factor(npair)
    if all(e == 1 for e in f.values()): sqfree_pairs += 1
    lam = abs(math.atan2(gmul(z, gconj(w))[1], gmul(z, gconj(w))[0]))
    liou = math.sqrt(2)/math.sqrt(npair)
    assert lam >= liou*(1-1e-12)
    # Matveev: logs log(pi_p/pibar_p) for p | npair, plus log i
    # Matveev's A_j = max(D h(alpha_j), |Log alpha_j|, 0.16) with D = 2, h(pi/pibar) = (log p)/2,
    # |Log(pi/pibar)| = 2|arg pi| <= pi; for log i: A = pi/2.  B is replaced by the lower bound
    # max(a_p, 4) (a larger B only weakens Matveev's bound, so this is conservative).
    A_prod = math.pi/2
    for p in f:
        g = gauss_prime(p)
        A_prod *= max(math.log(p), abs(2*math.atan2(g[1], g[0])), 0.16)
    B = max(list(f.values()) + [4])
    mv = matveev_log10(len(f)+1, A_prod, B)
    rows.append((npair, f, lam, liou, mv))
rows.sort()
log(f"  45 pairs; pair norms with all exponents 1: {sqfree_pairs}")
log("  pair norm n | factorisation | angle | sqrt2/sqrt n | ratio | log10 Matveev bound")
for npair, f, lam, liou, mv in rows[:8] + rows[-3:]:
    log(f"  {npair:>12} {str(f):<34} {lam:.3e} {liou:.3e} {lam/liou:7.2f} {mv:12.1f}")
log(f"  min over pairs of angle/(sqrt2/sqrt n) = {min(r[2]/r[3] for r in rows):.3f};"
    f" max = {max(r[2]/r[3] for r in rows):.1f}")
log(f"  Matveev's bound is weaker than the elementary bound for all 45 pairs: "
    f"{all(r[4] < math.log10(r[3]) for r in rows)}")

# (E') optimality family u = A + (A-1)i, w = i*ubar: angle = 2 arcsin(1/sqrt(2n))
worstrel = 0.0
for A in range(2, 2000):
    u = (A, A-1); n = gnorm(u)
    z = u; w = gmul((0, 1), gconj(u))
    lam = abs(math.atan2(gmul(z, gconj(w))[1], gmul(z, gconj(w))[0]))
    worstrel = max(worstrel, abs(lam/(math.sqrt(2)/math.sqrt(n)) - 1))
log(f"(E') pairs (u, i*ubar), u = A+(A-1)i, 2<=A<2000: max |angle/(sqrt2 n^(-1/2)) - 1| = {worstrel:.2e}")

# (F) product clusters
log("(F) product clusters with 3 coprime blocks: C >= sqrt2 * N^(1/12)")
worst = 1e9; tested = 0
for trial in range(3000):
    ps = random.sample(split_primes, random.randint(3, 6))
    random.shuffle(ps)
    cuts = sorted(random.sample(range(1, len(ps)), 2))
    blocks = [ps[:cuts[0]], ps[cuts[0]:cuts[1]], ps[cuts[1]:]]
    ab = []
    for bl in blocks:
        # two points of the block circle: product with random orientations, differing somewhere
        while True:
            e1 = [random.randint(0, 1) for _ in bl]; e2 = [random.randint(0, 1) for _ in bl]
            if e1 != e2: break
        def pt(e):
            z = (1, 0)
            for p, s in zip(bl, e):
                g = gauss_prime(p); z = gmul(z, g if s else gconj(g))
            return z
        ab.append((pt(e1), pt(e2)))
    pts = []
    for eps in itertools.product([0, 1], repeat=3):
        z = (1, 0)
        for (a, b), s in zip(ab, eps): z = gmul(z, b if s else a)
        pts.append(z)
    NN = gnorm(pts[0]); RR = math.sqrt(NN)
    ang = sorted(math.atan2(z[1], z[0]) % (2*math.pi) for z in pts)
    gaps = [(ang[(i+1) % 8] - ang[i]) % (2*math.pi) for i in range(8)]
    arc = RR*(2*math.pi - max(gaps))
    Cval = arc/math.sqrt(RR)
    ratio = Cval/(math.sqrt(2)*NN**(1/12))
    worst = min(worst, ratio); tested += 1
    assert ratio >= 1, (blocks, Cval)
log(f"  {tested} random block triples: min C/(sqrt2 N^(1/12)) = {worst:.3f} (lemma requires >= 1)")

# (G) rational directions
log("(G) arcs near Gaussian directions w: count <= 2|w|(kappa+C)^2+2")
viol = 0; checked = 0
Ns = []
for trial in range(400):
    ps = random.sample(split_primes[:25], random.randint(3, 5))
    NN = 1
    for p in ps: NN *= p
    Ns.append(NN)
for NN in Ns:
    if NN > 4*10**8: continue
    pts = circle_points(NN); RR = math.sqrt(NN)
    for w in [(1, 0), (1, 1), (2, 1), (1, 2), (3, 1), (3, 2), (2, 3), (1, 3)]:
        aw = math.atan2(w[1], w[0])
        for kappa, C in [(1.0, 1.0), (2.0, 1.0), (1.0, 3.0)]:
            # points whose angle is within (kappa + C)/sqrt(R) of arg w: superset of any qualifying arc
            near = [z for z in pts if abs(((math.atan2(z[1], z[0]) - aw + math.pi) % (2*math.pi)) - math.pi)
                    <= (kappa + C)/math.sqrt(RR)]
            bound = 2*math.sqrt(gnorm(w))*(kappa + C)**2 + 2
            checked += 1
            if len(near) > bound: viol += 1
log(f"  {checked} (circle, w, kappa, C) cases; violations of the bound: {viol}")

with open(__file__.replace('.py', '_output.txt'), 'w') as fh:
    fh.write('\n'.join(out) + '\n')
