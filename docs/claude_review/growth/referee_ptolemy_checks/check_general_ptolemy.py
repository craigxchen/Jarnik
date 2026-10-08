"""Exact verification of the general Ptolemy matching relation (Part 1 of ptolemy.md).

For actual lattice points z on x^2+y^2=N (N odd, all prime factors split, arbitrary
exponents, arbitrary Gaussian units) we check, for every quadruple in cyclic order:

  (G1) z_i = eps_i g_ij h_ij,  z_j = eps_j g_ij conj(h_ij)   (pair factorization)
  (G2) z_i - z_j = g_ij v_ij,  v_ij = eps_i h_ij - eps_j conj(h_ij)
  (G3) g12 g34 = Gamma X1,  g13 g24 = Gamma X2,  g14 g23 = Gamma X3  (exact, in Z[i])
       with X_k = prod_p p^(b3-b2) over primes whose within-side matching is M_k
  (G4) X2 v13 v24 = X1 v12 v34 + X3 v14 v23            (exact Gaussian identity)
  (G5) v12v34, v14v23, v13v24 are positive integer multiples a,b,c of ONE primitive omega,
       and c X2 = a X1 + b X3 (positive integers); X1,X2,X3 pairwise coprime
  (G6) threshold-layer description: X_k = prod over layers (p,l) whose cut restricted to
       the quadruple is 2|2 of type M_k
  (G7) common-unit class: |v_ij| = 2|Im h_ij| and the anchored half-angle residues
       t_ij = Im(conj(P_i)P_j)/G_ij, Y_i = Im(P_i) satisfy |t_ij| = |Im h_ij|,
       |Y_i| = |Im h_0i|, G_ij = Norm gcd(P_i,P_j) = prod_{T contains i,j} n_T,
       and |t_ac t_bd| X_{ac|bd} = |t_ab t_cd| X_{ab|cd} + |t_ad t_bc| X_{ad|bc}.
All checks are exact integer / Gaussian-integer arithmetic; angles are used only to sort.
"""
import random, math, itertools, sys
from math import gcd
from gi import *

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)

PRIMES = split_primes(120)
PI = {p: gaussian_prime_above(p) for p in PRIMES}

def point(unit, alloc, primes, exps):
    z = UNITS[unit]
    for p, a, e in zip(primes, alloc, exps):
        z = mul(z, mul(gpow(PI[p], a), gpow(conj(PI[p]), e - a)))
    return z

def pair_data(Ai, Aj, primes, exps):
    g = (1, 0); h = (1, 0)
    for p, ai, aj, e in zip(primes, Ai, Aj, exps):
        g = mul(g, mul(gpow(PI[p], min(ai, aj)), gpow(conj(PI[p]), e - max(ai, aj))))
        if ai > aj:
            h = mul(h, gpow(PI[p], ai - aj))
        elif aj > ai:
            h = mul(h, gpow(conj(PI[p]), aj - ai))
    return g, h

MATCHINGS = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]  # M1, M2, M3 on positions

def matching_X(allocs, primes, exps):
    """X_k from sorted valuations; also via threshold layers (G6)."""
    X = [1, 1, 1]
    Xl = [1, 1, 1]
    Gam = (1, 0)
    for idx, (p, e) in enumerate(zip(primes, exps)):
        vals = [allocs[q][idx] for q in range(4)]
        order = sorted(range(4), key=lambda q: vals[q])
        b = [vals[q] for q in order]
        Gam = mul(Gam, mul(gpow(PI[p], b[0] + b[1]), gpow(conj(PI[p]), 2*e - b[2] - b[3])))
        m = b[2] - b[1]
        if m > 0:
            low = set(order[:2])
            for k, (P1, P2) in enumerate(MATCHINGS):
                if set(P1) == low or set(P2) == low:
                    X[k] *= p**m
        # threshold layers: layer l (1..e) has side S_l = {q: a_q >= l}
        for l in range(1, e + 1):
            S = frozenset(q for q in range(4) if vals[q] >= l)
            if len(S) == 2:
                for k, (P1, P2) in enumerate(MATCHINGS):
                    if S == frozenset(P1) or S == frozenset(P2):
                        Xl[k] *= p
    return X, Xl, Gam

def check_quadruple(pts, primes, exps, stats):
    # pts: list of (unit, alloc, z) in cyclic (angular) order
    units = [UNITS[u] for u, _, _ in pts]
    allocs = [a for _, a, _ in pts]
    zs = [z for _, _, z in pts]
    g = {}; v = {}; h = {}
    for i, j in itertools.combinations(range(4), 2):
        gij, hij = pair_data(allocs[i], allocs[j], primes, exps)
        assert mul(units[i], mul(gij, hij)) == zs[i]                      # G1
        assert mul(units[j], mul(gij, conj(hij))) == zs[j]
        vij = sub(mul(units[i], hij), mul(units[j], conj(hij)))
        assert mul(gij, vij) == sub(zs[i], zs[j])                         # G2
        g[i, j] = gij; v[i, j] = vij; h[i, j] = hij
    X, Xl, Gam = matching_X(allocs, primes, exps)
    assert X == Xl                                                         # G6
    assert mul(g[0, 1], g[2, 3]) == smul(X[0], Gam)                        # G3
    assert mul(g[0, 2], g[1, 3]) == smul(X[1], Gam)
    assert mul(g[0, 3], g[1, 2]) == smul(X[2], Gam)
    alpha = mul(v[0, 1], v[2, 3]); beta = mul(v[0, 3], v[1, 2]); gam = mul(v[0, 2], v[1, 3])
    assert smul(X[1], gam) == add(smul(X[0], alpha), smul(X[2], beta))     # G4
    om = primitive_direction(gam)
    for w in (alpha, beta):
        assert primitive_direction(w) == om                               # same ray (positive multiples)
    a = content(alpha); b = content(beta); c = content(gam)
    assert c * X[1] == a * X[0] + b * X[2] and a > 0 and b > 0 and c > 0  # G5
    assert gcd(X[0], X[1]) == 1 and gcd(X[0], X[2]) == 1 and gcd(X[1], X[2]) == 1
    stats['quads'] += 1
    if len(set(units)) > 1:
        stats['mixed'] += 1
    if max(X) > 1:
        stats['nontrivialX'] += 1
    if any(e > 1 for e in exps):
        stats['prime_powers'] += 1
    return X, (a, b, c), om, h

def anchored_check(cluster, primes, exps, stats):
    """(G7) for a common-unit cluster, anchor = cluster[0]."""
    allocs = [a for _, a, _ in cluster]
    k = len(cluster)
    m = k - 1
    a0 = allocs[0]
    P = [(1, 0)]
    for i in range(1, k):
        Pi = (1, 0)
        for idx, p in enumerate(primes):
            d = allocs[i][idx] - a0[idx]
            Pi = mul(Pi, gpow(PI[p], d) if d > 0 else gpow(conj(PI[p]), -d))
        P.append(Pi)
    # blocks n_T from threshold layers relative to the anchor
    nT = {}
    for idx, (p, e) in enumerate(zip(primes, exps)):
        for l in range(1, e + 1):
            f0 = a0[idx] >= l
            T = frozenset(i for i in range(1, k) if (allocs[i][idx] >= l) != f0)
            if T:
                nT[T] = nT.get(T, 1) * p
    def Gprod(S):
        r = 1
        for T, n in nT.items():
            if S <= T:
                r *= n
        return r
    minor = {}
    resid = {}
    for i in range(k):
        for j in range(i + 1, k):
            mij = (conj(P[i])[0]*P[j][1] + conj(P[i])[1]*P[j][0])   # Im(conj(P_i) P_j)
            minor[i, j] = mij
            if i == 0:
                resid[i, j] = P[j][1]
                assert mij == P[j][1]
            else:
                Gij = norm(ggcd(P[i], P[j]))
                assert Gij == Gprod(frozenset([i, j]))                 # G_ij = prod_{T>=ij} n_T
                assert mij % Gij == 0
                resid[i, j] = mij // Gij
            # |Im h_ij| equals |residual|
            _, hij = pair_data(allocs[i], allocs[j], primes, exps)
            assert abs(hij[1]) == abs(resid[i, j]) and resid[i, j] != 0
            gij, _ = pair_data(allocs[i], allocs[j], primes, exps)
            vij = sub(mul(UNITS[cluster[i][0]], hij), mul(UNITS[cluster[j][0]], conj(hij)))
            assert norm(vij) == 4 * hij[1]**2
    # N(P_i) = prod_{T contains i} n_T
    for i in range(1, k):
        assert norm(P[i]) == Gprod(frozenset([i]))
    # anchored ordered matching identity for every 4-subset in angular order
    ang = [math.atan2(z[1], z[0]) for _, _, z in cluster]
    for Q in itertools.combinations(range(k), 4):
        Qs = sorted(Q, key=lambda q: -ang[q])
        if max(ang[q] for q in Q) - min(ang[q] for q in Q) > math.pi:
            continue
        sub_alloc = [allocs[q] for q in Qs]
        X, Xl, _ = matching_X(sub_alloc, primes, exps)
        r = lambda a, b: abs(resid[min(a, b), max(a, b)])
        A, B, C, D = Qs
        assert r(A, C)*r(B, D)*X[1] == r(A, B)*r(C, D)*X[0] + r(A, D)*r(B, C)*X[2]
        stats['anchored_quads'] += 1

def random_cluster(k, r, maxe, mixed):
    primes = random.sample(PRIMES, r)
    exps = [random.randint(1, maxe) for _ in primes]
    seen = set(); cl = []
    tries = 0
    while len(cl) < k and tries < 1000:
        tries += 1
        u = random.randrange(4) if mixed else 0
        alloc = tuple(random.randint(0, e) for e in exps)
        if (u, alloc) in seen:
            continue
        seen.add((u, alloc))
        cl.append((u, alloc, point(u, alloc, primes, exps)))
    return cl, primes, exps

def all_points_window(primes, exps, k):
    """all lattice points of the circle (unit class 0..3), sorted by angle; return shortest k-window"""
    pts = []
    for u in range(4):
        for alloc in itertools.product(*[range(e + 1) for e in exps]):
            pts.append((u, alloc, point(u, alloc, primes, exps)))
    pts.sort(key=lambda t: math.atan2(t[2][1], t[2][0]))
    best = None
    n = len(pts)
    for s in range(n):
        win = [pts[(s + q) % n] for q in range(k)]
        a0 = math.atan2(win[0][2][1], win[0][2][0]); a1 = math.atan2(win[-1][2][1], win[-1][2][0])
        span = (a1 - a0) % (2 * math.pi)
        if best is None or span < best[0]:
            best = (span, win)
    return best

stats = dict(quads=0, mixed=0, nontrivialX=0, prime_powers=0, anchored_quads=0, windows=0)
for trial in range(400):
    k = random.randint(4, 7)
    cl, primes, exps = random_cluster(k, random.randint(2, 5), 4, mixed=True)
    cl.sort(key=lambda t: -math.atan2(t[2][1], t[2][0]))
    for Q in itertools.combinations(range(len(cl)), 4):
        check_quadruple([cl[q] for q in Q], primes, exps, stats)
for trial in range(300):
    k = random.randint(4, 7)
    cl, primes, exps = random_cluster(k, random.randint(2, 5), 4, mixed=False)
    anchored_check(cl, primes, exps, stats)
# genuinely short windows on actual circles
for trial in range(40):
    r = random.randint(3, 5)
    primes = random.sample(PRIMES[:12], r)
    exps = [random.randint(1, 3) for _ in primes]
    span, win = all_points_window(primes, exps, 6)
    R = math.sqrt(math.prod(p**e for p, e in zip(primes, exps)))
    win = sorted(win, key=lambda t: -math.atan2(t[2][1], t[2][0]))
    for Q in itertools.combinations(range(6), 4):
        check_quadruple([win[q] for q in Q], primes, exps, stats)
    stats['windows'] += 1
print("all exact checks passed:", stats)
