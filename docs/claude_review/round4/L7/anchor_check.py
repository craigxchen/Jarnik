"""EXACT check of Theorem 1 (residue-free anchoring).

Random 7-point configurations in P^1(Q) with three points normalized to inf=(1,0), 0=(0,1), 1=(1,1).
Anchored rows by ptolemy Lemma 5.2: D = lcm|y_j|, P_j = (D x_j/|y_j|, sgn y_j) (primitive).
At every prime p dividing some det, compute the p-adic tree edge lengths l_T(p) from the Gromov
products g_ij = v_p(det(w_i,w_j)) of the ORIGINAL primitive vectors (min formula of ptolemy B(iv)),
and verify  v_p(det(P_i,P_j)) = sum_{T contains i,j} l_T(p)   (i.e. delta_p = 0)  and  Y_j = +-1.
Hence det(P_i,P_j) = +- prod_{T contains i,j} N_T exactly: H_Sigma = 0 for every configuration whose
trees are admissible (at most two nested chains).  Also tallies how often trees are non-admissible.
"""
import random, itertools
from math import gcd

def vp(n, p):
    n = abs(n); c = 0
    while n % p == 0: n //= p; c += 1
    return c

def primes_of(n):
    n = abs(n); ps = set(); d = 2
    while d*d <= n:
        while n % d == 0: ps.add(d); n //= d
        d += 1
    if n > 1: ps.add(n)
    return ps

def det(u, v): return u[0]*v[1] - u[1]*v[0]

def prim(u):
    g = gcd(u[0], u[1]); return (u[0]//g, u[1]//g)

def lcm(a, b): return a*b//gcd(a, b)

m = 6
rows = list(range(1, m+1))
Ts = [frozenset(c) for r in range(2, m) for c in itertools.combinations(rows, r)]  # visible: 2..m-1
TOP = frozenset(rows)

def edges(g, idx):
    """g: dict (i,j)->Gromov product, indices idx including 0.  Returns l_T for T subset rows, 2<=|T|<=m (incl. top)."""
    l = {}
    for T in Ts + [TOP]:
        out = [x for x in idx if x not in T]
        best = None
        for a, b in itertools.combinations(sorted(T), 2):
            for c, d in itertools.combinations(sorted(out), 2):
                G = lambda x, y: g[(min(x, y), max(x, y))]
                v = min(G(a, b) + G(c, d) - G(a, c) - G(b, d), G(a, b) + G(c, d) - G(a, d) - G(b, c))
                best = v if best is None else min(best, v)
        l[T] = max(0, best) if best is not None else 0
    return l

def chains_ok(Tset):
    """at most two nested chains on disjoint row sets"""
    Tl = sorted(Tset, key=len)
    chains = []
    for T in Tl:
        placed = False
        for ch in chains:
            if ch[-1] <= T:
                ch.append(T); placed = True; break
        if not placed: chains.append([T])
    if len(chains) > 2: return False
    if len(chains) == 2 and (chains[0][-1] & chains[1][-1]): return False
    return True

random.seed(1)
nconf = 0; nprime = 0; bad = 0; nonadm = 0; topnz = 0
for trial in range(400):
    pts = [(1, 0), (0, 1), (1, 1)]
    while len(pts) < 7:
        q = prim((random.randint(-60, 60), random.randint(1, 60)))
        if all(det(q, r) != 0 for r in pts): pts.append(q)
    w = pts  # w_0 = (1,0) is the anchor; rows 1..6
    y = [det(w[0], w[j]) for j in range(7)]
    D = 1
    for j in range(1, 7): D = lcm(D, abs(y[j]))
    P = [(1, 0)] + [prim((D*w[j][0]//abs(y[j]), (1 if y[j] > 0 else -1))) for j in range(1, 7)]
    # Y_j
    assert all(abs(det(P[0], P[j])) == 1 for j in range(1, 7))
    allp = set()
    for i, j in itertools.combinations(range(7), 2): allp |= primes_of(det(w[i], w[j]))
    for p in allp:
        g = {(i, j): vp(det(w[i], w[j]), p) for i, j in itertools.combinations(range(7), 2)}
        l = edges(g, list(range(7)))
        for i, j in itertools.combinations(rows, 2):
            lhs = vp(det(P[i], P[j]), p)
            rhs = sum(l[T] for T in Ts + [TOP] if i in T and j in T)
            if lhs != rhs: bad += 1
        if l[TOP] > 0: topnz += 1
        pos = [T for T in Ts if l[T] > 0]
        if not chains_ok(pos): nonadm += 1
        nprime += 1
    nconf += 1
print('configurations', nconf, 'prime checks', nprime, 'mismatches', bad,
      'primes with top-block length > 0:', topnz, 'non-admissible trees:', nonadm)
