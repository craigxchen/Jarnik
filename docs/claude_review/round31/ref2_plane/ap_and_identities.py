"""Referee ref2_plane: independent checks.
(a) AP direction enumeration: no 5 real points whose 10 segment slopes are exactly {0,...,9}.
(b) primality of 34913+2100r, Gaussian blocks, cofactors t_ij, identities (II.10),(II.6),(8), Pluecker forms.
"""
from fractions import Fraction as Fr
from itertools import permutations, combinations
import random

# (a) Independent enumeration: brute force over all ways, not only the anchor normalisation.
# Points A_1..A_5; the unique slope-0 segment normalised to (0,0),(1,0).
cnt = 0; sols = 0
S = list(range(1, 10))
seen = set()
for six in permutations(S, 6):
    pts = []
    for k in range(3):
        a, b = six[2*k], six[2*k+1]
        # point P with slope a from (0,0), slope b from (1,0)
        x = Fr(b, b - a); y = a * x
        pts.append((x, y))
    key = frozenset(pts)
    if key in seen:
        continue
    seen.add(key)
    cnt += 1
    rest = set(S) - set(six)
    internal = []
    ok = True
    for P, Q in combinations(pts, 2):
        if P[0] == Q[0]:
            ok = False; break
        internal.append((Q[1]-P[1])/(Q[0]-P[0]))
    if ok and sorted(internal) == sorted(Fr(r) for r in rest):
        sols += 1
print("AP enumeration: distinct triples", cnt, "solutions", sols)

# Second independent check: search over all real configurations via a different parametrisation:
# label edges by slopes; any 5-point config with distinct slopes {0..9}: pick the pair with slope 0 and slope 1 etc.
# Brute force over assignments of slopes to the 10 edges of K5 (10! / symmetries), solve linear system
# for the points; check consistency.  Use: A1=(0,0); A2 on slope s12 line at x=1; then others determined by
# lines from A1, A2; check remaining slopes.
from itertools import permutations as perms
def config(assign):
    # assign: dict edge->slope
    A = {1: (Fr(0), Fr(0)), 2: (Fr(1), Fr(assign[(1,2)]))}
    for k in (3,4,5):
        a = assign[(1,k)]; b = assign[(2,k)]
        if a == b: return None
        # y = a x ; y - y2 = b (x - 1)
        x = (A[2][1] - b) / (a - b); y = a*x
        A[k] = (x, y)
    for (i,j) in [(3,4),(3,5),(4,5)]:
        if A[i][0] == A[j][0]: return None
        if (A[j][1]-A[i][1])/(A[j][0]-A[i][0]) != assign[(i,j)]: return None
    return A
edges = [(i,j) for i in range(1,6) for j in range(i+1,6)]
found = 0; tested = 0
# fix slope of edge (1,2) = 0 wlog? Not wlog under labels; enumerate all 10! assignments of slopes to edges
# but use symmetry: relabel so that edge with slope 0 is (1,2): then 8!*... we just do full 10!/ (cheap enough?) 3.6M
for perm in perms(range(10)):
    if perm[0] != 0:  # edge (1,2) gets slope 0 (wlog by relabelling marks)
        break
    assign = dict(zip(edges, perm))
    tested += 1
    if config(assign) is not None:
        found += 1
print("full assignment search (slope 0 on edge 12):", tested, "assignments, configs found:", found)

# (b) primes 34913+2100r
def isprime(n):
    if n < 2: return False
    d = 2
    while d*d <= n:
        if n % d == 0: return False
        d += 1
    return True
ps = [34913 + 2100*r for r in range(10)]
print("primes:", ps, all(isprime(p) for p in ps), [p % 4 for p in ps])

def gauss_rep(p):
    for a in range(1, int(p**0.5)+1):
        b2 = p - a*a
        b = int(round(b2**0.5))
        if b*b == b2: return (a, b)
def gmul(x, y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def conj(x): return (x[0], -x[1])
def det(x, y): return x[0]*y[1]-x[1]*y[0]
def dot(x, y): return x[0]*y[0]+x[1]*y[1]

def check_tuple(norm_of_edge, rng):
    beta = {}
    for e, p in norm_of_edge.items():
        a, b = gauss_rep(p)
        g = rng.choice([(a, b), (b, a), (a, -b), (-a, b)])
        beta[e] = g
    P = {}
    for i in range(1, 6):
        v = (1, 0)
        for j in range(1, 6):
            if j != i:
                v = gmul(v, beta[tuple(sorted((i, j)))])
        P[i] = v
    n = {e: norm_of_edge[e] for e in norm_of_edge}
    t = {}; s = {}
    for i in range(1, 6):
        for j in range(1, 6):
            if i == j: continue
            e = tuple(sorted((i, j)))
            d = det(P[i], P[j]); dd = dot(P[i], P[j])
            assert d % n[e] == 0 and dd % n[e] == 0
            t[(i, j)] = d // n[e]; s[(i, j)] = dd // n[e]
    # II.10
    for i in range(1,6):
        for j in range(1,6):
            if i==j: continue
            prod = 1
            for a in range(1,6):
                if a not in (i,j):
                    prod *= n[tuple(sorted((i,a)))]*n[tuple(sorted((j,a)))]
            assert s[(i,j)]**2 + t[(i,j)]**2 == prod
    # II.6 for all ordered triples
    N = lambda a,b: n[tuple(sorted((a,b)))]
    for i,j,k in permutations(range(1,6),3):
        l,m = [x for x in range(1,6) if x not in (i,j,k)]
        A = t[(j,k)]**2*N(j,k)*N(i,l)*N(i,m) + t[(i,k)]**2*N(i,k)*N(j,l)*N(j,m) - t[(i,j)]**2*N(i,j)*N(k,l)*N(k,m)
        assert 2*t[(i,k)]*t[(j,k)]*s[(i,j)] == A
    # (8)
    tot = [Fr(0), Fr(0)]
    for i in range(1,6):
        pr = 1
        for j in range(1,6):
            if j != i: pr *= t[(i,j)]
        tot[0] += Fr(P[i][0], pr); tot[1] += Fr(P[i][1], pr)
    assert tot == [0, 0]
    # Pluecker on each quartet with M_ij = t_ij n_ij
    M = lambda a,b: t[(a,b)]*N(a,b)
    for (i,j,k,l) in combinations(range(1,6),4):
        assert M(i,j)*M(k,l) - M(i,k)*M(j,l) + M(i,l)*M(j,k) == 0
    return t
rng = random.Random(1)
edge_rank = {e: r for r, e in enumerate(edges)}
for trial in range(5):
    perm = list(range(10)); rng.shuffle(perm)
    noe = {e: ps[perm[r]] for r, e in enumerate(edges)}
    t = check_tuple(noe, rng)
    T = max(abs(v) for v in t.values())
    print("AP tuple trial", trial, "all t nonzero:", all(v != 0 for v in t.values()), "T =", T,
          "sqrt(max n)=%.1f" % (max(ps)**0.5))
# random genuine tuples with random primes 1 mod 4: identities
primes1 = [p for p in range(5, 4000) if isprime(p) and p % 4 == 1]
for trial in range(200):
    pr = rng.sample(primes1, 10)
    check_tuple(dict(zip(edges, pr)), rng)
print("identities (II.10),(II.6),(8),Pluecker verified on 200 random genuine ten-pair tuples")
