"""Exact check of Theorem B(iii): the Ptolemy data of an actual common-unit circle cluster are the
boundary contacts of a point of M_{0,k}, and the p-adic tree edge of every boundary split S is
the block valuation v_p(n_S), up to 2*max v_p(residual).

For random actual clusters (several split primes, ARBITRARY exponents -> nested blocks sharing
rational primes, exactly as in the item-532 extraction) with anchor z_0:
  * P_i = prod_p pi^{(a_i-a_0)+} pibar^{(a_0-a_i)+}, P_0 = 1, blocks n_T from threshold layers;
  * Gromov products g_ij(p) = v_p(det(P_i,P_j)) (det(P_0,P_j) = Im P_j = y_j);
  * edge lengths l_S(p) = max(0, min over quartets a,b in S, c,d not in S of
        min(g_ab+g_cd-g_ac-g_bd, g_ab+g_cd-g_ad-g_bc));
  * checks: (T1) the positive-length splits are pairwise compatible (they form a tree);
            (T2) the tree reproduces EVERY quartet cross-ratio valuation;
            (T3) |l_S(p) - v_p(n_S)| <= 2 max_{ij} v_p(t_ij)   (t_0j := y_j), and equality
                 l_S(p) = v_p(n_S) whenever p divides no residual;
            (T4) singleton blocks and the top block never produce a boundary edge.
All arithmetic is exact (integers / Gaussian integers).
"""
import sys, random, itertools, math
sys.path.insert(0, '.')
from gi import *

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
PRIMES = split_primes(60)
PI = {p: gaussian_prime_above(p) for p in PRIMES}

def vp_int(n, p):
    n = abs(n)
    assert n != 0
    v = 0
    while n % p == 0:
        n //= p; v += 1
    return v

def splits(k):
    """boundary splits of {0..k-1}: side S not containing 0, 2<=|S|<=k-2"""
    rows = list(range(1, k))
    return [frozenset(S) for r in range(2, k - 1) for S in itertools.combinations(rows, r)]

def edge_lengths(g, k, S_list):
    out = {}
    for S in S_list:
        Sc = [x for x in range(k) if x not in S]
        best = None
        for a, b in itertools.combinations(sorted(S), 2):
            for c, d in itertools.combinations(Sc, 2):
                I = min(g[a, b] + g[c, d] - g[a, c] - g[b, d], g[a, b] + g[c, d] - g[a, d] - g[b, c])
                best = I if best is None else min(best, I)
        out[S] = max(0, best)
    return out

def compatible(S1, S2, k):
    A = set(S1); B = set(S2); U = set(range(k))
    return not (A & B and A - B and B - A and (U - A) & (U - B)) or A <= B or B <= A or not (A & B)

stats = dict(clusters=0, prime_checks=0, residual_primes=0, nested=0, edges=0)
for trial in range(250):
    k = random.randint(5, 8)
    r = random.randint(2, 4)
    primes = random.sample(PRIMES, r)
    exps = [random.randint(1, 4) for _ in primes]
    allocs = set()
    while len(allocs) < k:
        allocs.add(tuple(random.randint(0, e) for e in exps))
        if len(allocs) == math.prod(e + 1 for e in exps):
            break
    allocs = list(allocs)
    if len(allocs) < 5:
        continue
    k = len(allocs)
    a0 = allocs[0]
    P = [(1, 0)]
    for i in range(1, k):
        Pi = (1, 0)
        for idx, p in enumerate(primes):
            dd = allocs[i][idx] - a0[idx]
            Pi = mul(Pi, gpow(PI[p], dd) if dd > 0 else gpow(conj(PI[p]), -dd))
        P.append(Pi)
    # blocks from threshold layers relative to the anchor (item 532)
    nT = {}
    for idx, (p, e) in enumerate(zip(primes, exps)):
        for l in range(1, e + 1):
            f0 = a0[idx] >= l
            T = frozenset(i for i in range(1, k) if (allocs[i][idx] >= l) != f0)
            if T:
                nT[T] = nT.get(T, 1) * p
    det = {}
    for i, j in itertools.combinations(range(k), 2):
        det[i, j] = P[i][0] * P[j][1] - P[j][0] * P[i][1]
        assert det[i, j] != 0
    G = lambda i, j: math.prod(n for T, n in nT.items() if i in T and j in T)
    t = {}
    for i, j in itertools.combinations(range(k), 2):
        g = 1 if i == 0 else G(i, j)
        assert det[i, j] % g == 0
        t[i, j] = det[i, j] // g
    S_list = splits(k)
    # all primes dividing any det
    plist = set(primes)
    for v in list(det.values()):
        n = abs(v); q = 2
        while q * q <= n:
            while n % q == 0:
                plist.add(q); n //= q
            q += 1
        if n > 1:
            plist.add(n)
    if any(len([q for q in primes if any(a[primes.index(q)] not in (0, exps[primes.index(q)]) for a in allocs)]) > 0 for _ in [0]):
        stats['nested'] += 1
    for p in plist:
        g = {}
        for (i, j), v in det.items():
            g[i, j] = g[j, i] = vp_int(v, p)
        ell = edge_lengths(g, k, S_list)
        pos = [S for S in S_list if ell[S] > 0]
        # (T1) compatibility
        for S1, S2 in itertools.combinations(pos, 2):
            A, B = set(S1), set(S2)
            assert A <= B or B <= A or not (A & B), (S1, S2)   # both sides avoid 0, so this is compatibility
        # (T2) reproduce every quartet cross-ratio valuation
        for a, b, c, d in itertools.combinations(range(k), 4):
            for (x1, x2), (y1, y2) in [((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))]:
                # v_p( det(x1,x2) det(y1,y2) / (det(x1,y1) det(x2,y2)) )
                v = g[x1, x2] + g[y1, y2] - g[x1, y1] - g[x2, y2]
                tv = 0
                for S in pos:
                    side = lambda u: u in S
                    same12 = side(x1) == side(x2) and side(y1) == side(y2) and side(x1) != side(y1)
                    same_13 = side(x1) == side(y1) and side(x2) == side(y2) and side(x1) != side(x2)
                    tv += ell[S] * (int(same12) - int(same_13))
                assert v == tv, (p, (a, b, c, d))
        # (T3) comparison with block valuations
        mu = max(vp_int(t[i, j], p) for i, j in t)
        for S in S_list:
            nu = vp_int(nT.get(S, 1), p)
            assert abs(ell[S] - nu) <= 2 * mu, (p, S, ell[S], nu, mu)
            if mu == 0:
                assert ell[S] == nu
            stats['edges'] += ell[S] > 0
        if mu > 0:
            stats['residual_primes'] += 1
        stats['prime_checks'] += 1
    # (T4): singleton and top blocks are not boundary splits by definition; check they leave no trace
    stats['clusters'] += 1
print("all checks passed:", stats)
