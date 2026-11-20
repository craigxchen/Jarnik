"""Independent referee check of the Paley squarefree relaxation (Astra IX.4).
Exact integer arithmetic throughout. Own prime choice (see choose_primes)."""
import sys, math
from itertools import combinations
from math import gcd, prod, isqrt
from sympy import primerange, Matrix

def legendre(a, m):
    a %= m
    if a == 0: return 0
    return 1 if pow(a, (m-1)//2, m) == 1 else -1

def cut_checks(m):
    K = [[legendre(v-i, m) for v in range(m)] for i in range(m)]
    S = [[K[i][v] - (1 if i == v else 0) for v in range(m)] for i in range(m)]
    SM = Matrix(S)
    G1 = SM.T*SM; G2 = SM*SM.T
    target = (m+1)*Matrix.eye(m) - Matrix.ones(m, m)
    assert G1 == target and G2 == target
    one = Matrix.ones(m, 1)
    assert SM*one == -one and SM.T*one == -one
    assert all(K[i][v] == -K[v][i] for i in range(m) for v in range(m))
    agree = {sum(S[i][v] == S[j][v] for v in range(m)) for i, j in combinations(range(m), 2)}
    A = [[(1+S[i][v])//2 for v in range(m)] for i in range(m)]
    wts = {sum(r) for r in A}
    ovl = {sum(A[i][v]*A[j][v] for v in range(m)) for i, j in combinations(range(m), 2)}
    AAt = (Matrix(A)*Matrix(A).T).applyfunc(lambda x: x % 2)
    rows = [sum(b << v for v, b in enumerate(r)) for r in A]
    xors = {rows[i] ^ rows[j] for i, j in combinations(range(m), 2)}
    # column sums of S: plus/minus counts, agreeing pairs per column
    colagree = {sum(S[i][v] == S[j][v] for i, j in combinations(range(m), 2)) for v in range(m)}
    print(f"m={m} m%8={m%8}: S^tS=SS^t=(m+1)I-J ok; S1=S^t1=-1 ok; K skew ok; "
          f"row agreements {agree}; A row weights {wts}; overlaps {ovl}; "
          f"AA^t=I mod 2: {AAt == Matrix.eye(m)}; distinct pair xors {len(xors)} of {m*(m-1)//2}; "
          f"agreeing pairs per column {colagree} (min possible {((m-1)//2)**2})")
    return S, A, rows

def choose_primes(m, P):
    """Primes 1 mod 4 in [P,2P]; products in bins of ratio exp(1/J); greedy matching in best bin."""
    J = math.ceil(m*math.log(m))
    pr = [p for p in primerange(P, 2*P+1) if p % 4 == 1]
    # bins by floor(J*log(product/P^2))
    from collections import defaultdict
    bins = defaultdict(list)
    for a, b in combinations(pr, 2):
        bins[math.floor(J*math.log(a*b/(P*P)))].append((a, b))
    best = None
    for key, edges in bins.items():
        deg = defaultdict(int)
        for a, b in edges: deg[a] += 1; deg[b] += 1
        used = set(); match = []
        for a, b in edges:  # greedy
            if a not in used and b not in used:
                match.append((a, b)); used.update((a, b))
        if best is None or len(match) > len(best[1]):
            best = (key, match, len(edges), max(deg.values()))
    key, match, ne, md = best
    return J, len(pr), key, match, ne, md

def run(m, P):
    S, A, rows = cut_checks(m)
    J, npr, key, match, ne, md = choose_primes(m, P)
    print(f"P={P}, c0=P/(m log^2 m)={P/(m*math.log(m)**2):.2f}, J={J}, P/J={P/J:.1f}, primes 1 mod 4 in [P,2P]: {npr}; "
          f"best bin #{key}: {ne} products, max degree {md}, greedy matching {len(match)}; "
          f"edges/(2*maxdeg)={ne/(2*md):.1f}; P/log^2P={P/math.log(P)**2:.1f}")
    if len(match) < m:
        print("  not enough"); return None
    match = sorted(match, key=lambda e: e[0]*e[1])[:m]
    Q = [a*b for a, b in match]
    primes = sorted(p for e in match for p in e)
    assert len(set(primes)) == 2*m
    T = min(Q)
    assert max(Q) < T*math.exp(1/J) * (1+1e-12)
    q = prod(Q); W = math.log(q)
    k = P//(4*m)
    w = [k*i + (1 if i == m-1 else 0) for i in range(m)]
    d = {}; h = {}
    for i, j in combinations(range(m), 2):
        H = prod(Q[v] for v in range(m) if S[i][v] == S[j][v])
        t = w[j]-w[i]
        assert 0 < abs(t) < min(primes)
        d[i, j] = H*t; h[i, j] = H
        assert gcd(d[i, j], q) == H
    r = isqrt(q)
    maxd = max(abs(x) for x in d.values())
    assert maxd*maxd < q
    ratio_h = [hv*hv*T/q for hv in h.values()]  # (h/sqrt(q/T))^2
    print(f"  s={2*m}, log q={W:.1f}, 2 m log m={2*m*math.log(m):.1f}, log q/loglog q={W/math.log(W):.2f} (=> ratio m/that {m/(W/math.log(W)):.3f}); "
          f"max|d|/sqrt q={maxd/r:.4f}; (h/sqrt(q/T)) in [{math.sqrt(min(ratio_h)):.5f},{math.sqrt(max(ratio_h)):.5f}]; 1/log m={1/math.log(m):.3f}")
    gt = 0
    for i, j in combinations(range(m), 2): gt = gcd(gt, w[j]-w[i])
    meanlogt = sum(math.log(abs(w[j]-w[i])) for i, j in combinations(range(m), 2))/(m*(m-1)/2)
    print(f"  joint gcd of t: {gt}; mean log|t|={meanlogt:.3f}; log m={math.log(m):.3f}; log k={math.log(k):.3f}")
    # Sidon signatures in F_2^{2m} (each column duplicated)
    sig = [sum(((A[i][v]) << (2*v)) | ((A[i][v]) << (2*v+1)) for v in range(m)) for i in range(m)]
    sums = {sig[i] ^ sig[j] for i, j in combinations(range(m), 2)}
    print(f"  distinct unordered pair sums of signatures in F_2^{2*m}: {len(sums)} (binom={m*(m-1)//2}); D={1+len(sums)}")
    # triangle identities, unordered triples i<j<l: d_ij + d_jl = d_il
    fails = sum(d[i, j] + d[j, l] != d[i, l] for i, j, l in combinations(range(m), 3))
    print(f"  failed signed triangle identities (unordered triples i<j<l): {fails} of {math.comb(m,3)}")
    # mod p (p | q) triangle identities, and the pair residue test: separated pairs must have d^2 = const mod p
    pmap = {}
    for v, (a, b) in enumerate(match[:m]): pmap[a] = v; pmap[b] = v
    modfail_min = None; distinct_sq = []
    for p, v in pmap.items():
        f = sum((d[i, j] + d[j, l] - d[i, l]) % p != 0 for i, j, l in combinations(range(m), 3))
        modfail_min = f if modfail_min is None else min(modfail_min, f)
        sq = {(d[i, j]*d[i, j]) % p for i, j in combinations(range(m), 2) if S[i][v] != S[j][v]}
        distinct_sq.append(len(sq))
    print(f"  mod p triangle failures, min over the {2*m} primes p|q: {modfail_min}; "
          f"distinct values of d^2 mod p over separated pairs: min {min(distinct_sq)}, max {max(distinct_sq)} (genuine twin: 1)")
    for l in [3, 5, 7]:
        f = sum((d[i, j] + d[j, ll] - d[i, ll]) % l != 0 for i, j, ll in combinations(range(m), 3))
        print(f"  mod {l} triangle failures: {f}")
    # local valuation witnesses at p | q: marks b + p w
    for p, v in list(pmap.items())[:4]:
        loc = [(1+S[i][v])//2 + p*w[i] for i in range(m)]
        def vp(x, p):
            c = 0
            while x % p == 0: x //= p; c += 1
            return c
        ok = all(vp(loc[j]-loc[i], p) == vp(d[i, j], p) for i, j in combinations(range(m), 2))
        assert ok
    # valuations at small l not dividing q: equal to those of w differences
    for l in primerange(2, 200):
        assert all((d[i, j] % l == 0) == ((w[j]-w[i]) % l == 0) for i, j in combinations(range(m), 2))
    print("  local valuation witnesses verified (p|q: b+pw, first 4 primes; l<200: w)")
    # nonrealizability margin: q/(Q_u Q_v) versus (m+1)*sqrt q
    margin = min(q//(Q[u]*Q[v]) for u, v in combinations(range(m), 2)) > (m+1)*r + 1
    # column identity a^t S_w = 0
    for u, v in [(0, 1), (2, 5)]:
        a = [S[i][u]-S[i][v] for i in range(m)]
        assert sum(a) == 0 and all(sum(a[i]*S[i][x] for i in range(m)) == 0 for x in range(m) if x not in (u, v))
        assert sum(abs(x) for x in a) == m+1
    print(f"  a^t S_w = 0 checked; sum|a| = m+1; q/(Q_uQ_v) > (m+1) sqrt q for all u,v: {margin}; "
          f"log(q/(Q_uQ_v))/log q = {(W-2*math.log(max(Q)))/W:.3f}")
    # Vandermonde/larger-sieve slack: binom(m,2) log sqrt q - sum log|d|  vs  log superfactorial
    lhs = math.comb(m, 2)*W/2 - sum(math.log(abs(x)) for x in d.values())
    sf = sum(math.lgamma(j+1) for j in range(1, m))
    print(f"  pair-valuation slack sum log(sqrt q/|d_ij|) = {lhs:.1f}; log prod_(j<m) j! = {sf:.1f}; "
          f"sum log|t| = {sum(math.log(abs(w[j]-w[i])) for i,j in combinations(range(m),2)):.1f}")
    return q

if __name__ == "__main__":
    for m in [11, 19, 43]:
        cut_checks(m)
    m = int(sys.argv[1]); P = int(sys.argv[2])
    run(m, P)
