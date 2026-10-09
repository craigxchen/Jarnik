"""EXACT independent check of L7 Theorem 1 (referee).

Independent of the Gromov-product min-formula used in anchor_check.py: the p-adic tree is computed
from the CLUSTER PICTURE with the anchor (label 0) at infinity.  With w_0 = (1,0) the points are
a_j = x_j / y_j in Q (j = 1..6).  For a proper subset T of [6] with |T| >= 2 that is a p-adic cluster
(T = the set of a_j in some disc), its edge length is

    l_T = min_{i,j in T} v(a_i - a_j)  -  max_{i in T, k in [6]\\T} v(a_i - a_k)      (if > 0),

and l_T = 0 for non-clusters.  (The disc of T and the next larger disc.)  Then we check

    v_p(det(P_i,P_j)) == sum_{T contains i,j, |T|<=5} l_T(p)          (Theorem 1, delta_p = 0)

for the anchored rows P_j = (D x_j/|y_j|, sgn y_j), D = lcm |y_j|.

Stress families: (A) random small configurations normalized at inf,0,1 (as in L7);
(B) configurations with deep nested p-adic clusters at p in {2,3,5,7} (large valuations, many
non-admissible trees); (C) the SAME point sets WITHOUT the median normalization (labels 1,2 not at
0,1), to confirm that delta_p = 0 genuinely needs it (expect delta_p > 0 sometimes).
"""
import random, itertools
from fractions import Fraction as Fr
from math import gcd

def vp_int(n, p):
    n = abs(n); c = 0
    if n == 0: return 10**9
    while n % p == 0: n //= p; c += 1
    return c

def vp(q, p):
    q = Fr(q)
    if q == 0: return 10**9
    return vp_int(q.numerator, p) - vp_int(q.denominator, p)

def primes_of(n):
    n = abs(n); ps = set(); d = 2
    while d * d <= n:
        while n % d == 0: ps.add(d); n //= d
        d += 1
    if n > 1: ps.add(n)
    return ps

def det(u, v): return u[0]*v[1] - u[1]*v[0]
def prim(u):
    g = gcd(u[0], u[1]); return (u[0]//g, u[1]//g)
def lcm(a, b): return a*b//gcd(a, b)

ROWS = list(range(1, 7))
TS = [frozenset(c) for r in range(2, 6) for c in itertools.combinations(ROWS, r)]

def cluster_lengths(a, p):
    """a: dict j -> Fraction (finite points, anchor at infinity)."""
    l = {}
    for T in TS:
        inner = min(vp(a[i] - a[j], p) for i, j in itertools.combinations(sorted(T), 2))
        outer = max(vp(a[i] - a[k], p) for i in T for k in ROWS if k not in T)
        l[T] = max(0, inner - outer)
    return l

def anchored(w):
    y = [det(w[0], w[j]) for j in range(7)]
    D = 1
    for j in range(1, 7): D = lcm(D, abs(y[j]))
    P = [(1, 0)] + [prim((D*w[j][0]//abs(y[j]), 1 if y[j] > 0 else -1)) for j in range(1, 7)]
    return P

def check(w, primes=None):
    """returns (#prime checks, #mismatches, list of (p, delta) with delta != 0)."""
    P = anchored(w)
    assert all(abs(det(P[0], P[j])) == 1 for j in range(1, 7))
    a = {j: Fr(w[j][0], w[j][1]) for j in ROWS}
    allp = set()
    for i, j in itertools.combinations(range(7), 2): allp |= primes_of(det(w[i], w[j]))
    if primes is not None: allp |= set(primes)
    mism = 0; deltas = []
    for p in sorted(allp):
        l = cluster_lengths(a, p)
        ds = set()
        for i, j in itertools.combinations(ROWS, 2):
            lhs = vp_int(det(P[i], P[j]), p)
            rhs = sum(l[T] for T in TS if i in T and j in T)
            ds.add(lhs - rhs)
        if ds != {0}:
            if len(ds) == 1 and min(ds) > 0: deltas.append((p, ds.pop()))
            else: mism += 1
    return len(allp), mism, deltas

def nonadmissible(w, p):
    a = {j: Fr(w[j][0], w[j][1]) for j in ROWS}
    l = cluster_lengths(a, p)
    pos = [T for T in TS if l[T] > 0]
    # exact test: can pos be covered by <= 2 chains whose top elements are disjoint?
    for r in range(0, len(pos) + 1):
        pass
    # brute force 2-colouring
    n = len(pos)
    for mask in range(1 << n):
        c1 = [pos[k] for k in range(n) if mask >> k & 1]
        c2 = [pos[k] for k in range(n) if not mask >> k & 1]
        ok = True
        for ch in (c1, c2):
            ch.sort(key=len)
            for x, y in zip(ch, ch[1:]):
                if not x <= y: ok = False; break
            if not ok: break
        if ok and c1 and c2 and (c1[-1] & c2[-1]): ok = False
        if ok: return False
    return True

random.seed(20261009)
tot = [0, 0, 0]; nonadm = 0; prime_checks = 0
# (A) random normalized
for _ in range(300):
    pts = [(1, 0), (0, 1), (1, 1)]
    while len(pts) < 7:
        q = prim((random.randint(-200, 200), random.randint(1, 200)))
        if all(det(q, r) != 0 for r in pts): pts.append(q)
    n, m, d = check(pts); tot[0] += n; tot[1] += m; tot[2] += len(d)
print('(A) random normalized: prime checks', tot[0], 'mismatches', tot[1], 'nonzero delta_p', tot[2])

# (B) deep nested clusters at small primes, normalized at inf,0,1
totB = [0, 0, 0]; nonadmB = 0; nB = 0
for _ in range(300):
    p = random.choice([2, 3, 5, 7])
    pts = [(1, 0), (0, 1), (1, 1)]
    while len(pts) < 7:
        base = random.choice([0, 1, random.randint(-5, 5)])
        x = Fr(base) + sum(Fr(random.randint(-3, 3) * p**e) for e in range(1, random.randint(2, 7)))
        x = x / random.choice([1, 1, p, p**2, 3])
        q = prim((x.numerator, x.denominator))
        if all(det(q, r) != 0 for r in pts): pts.append(q)
    n, m, d = check(pts, primes=[p]); totB[0] += n; totB[1] += m; totB[2] += len(d)
    nB += 1; nonadmB += nonadmissible(pts, p)
print('(B) deep clusters: prime checks', totB[0], 'mismatches', totB[1], 'nonzero delta_p', totB[2],
      '| non-admissible trees at the cluster prime:', nonadmB, 'of', nB)

# (C) same kind of data but labels 1,2 NOT at 0,1 (only anchor at infinity)
totC = [0, 0, 0]; ex = None
for _ in range(300):
    p = random.choice([2, 3, 5])
    pts = [(1, 0)]
    while len(pts) < 7:
        x = Fr(random.randint(-3, 3)) + sum(Fr(random.randint(-2, 2) * p**e) for e in range(1, 5))
        x = x * p**random.randint(0, 3)
        q = prim((x.numerator, x.denominator))
        if all(det(q, r) != 0 for r in pts): pts.append(q)
    n, m, d = check(pts, primes=[p]); totC[0] += n; totC[1] += m; totC[2] += len(d)
    if d and ex is None: ex = (pts, d)
print('(C) no median normalization: prime checks', totC[0], 'non-constant-offset mismatches', totC[1],
      'primes with delta_p > 0:', totC[2])
print('    example:', ex)
