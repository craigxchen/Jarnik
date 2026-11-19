"""All monic points on Astra's parent N (V.1, small d) and the largest c=0 cluster among them,
to show the family's points are a small part of the parent's point set.
usage: fullparent.py d L"""
import sys, itertools
import networkx as nx
from f4poly import *

d, L = int(sys.argv[1]), int(sys.argv[2])
Q = 1 << d
T = 1 << L
Y_factors = [linfac(T, v) for v in range(Q) if v not in (0, 1)]
N4 = prod(Y_factors)
assert all(x in (0, 1) for x in N4)
N = sum(1 << i for i, x in enumerate(N4) if x)   # as F_2[t] int
D = deg(N4) // 2

def f2divmod(a, b):
    q = 0
    db = f2deg(b)
    while a and f2deg(a) >= db:
        s = f2deg(a) - db
        q ^= 1 << s
        a ^= b << s
    return q, a

def f2factor(a):
    out = {}
    p = 2
    while f2deg(p) * 2 <= f2deg(a):
        while True:
            q, r = f2divmod(a, p)
            if r:
                break
            out[p] = out.get(p, 0) + 1
            a = q
        p += 1
    if a > 1:
        out[a] = out.get(a, 0) + 1
    return out  # (trial division by all polys, so factors found are irreducible)

fac = f2factor(N)
print(f"d={d} L={L}: deg N = {f2deg(N)}, D = {D}; F_2 factorisation degrees/exponents:",
      sorted((f2deg(p), e) for p, e in fac.items()))

def f4_split(P):
    e = f2deg(P)
    PP = from_f2(P)
    for coeffs in itertools.product(range(4), repeat=e // 2):
        f = tuple(coeffs) + (1,)
        q, r = divmod_(PP, f)
        if not r:
            return f, q
    raise ValueError

choices = []
for P, k in fac.items():
    if f2deg(P) % 2:  # inert
        assert k % 2 == 0
        choices.append([prod([from_f2(P)] * (k // 2))])
    else:
        f, fb = f4_split(P)
        assert fb == conj(f)
        choices.append([mul(prod([f] * i), prod([fb] * (k - i))) for i in range(k + 1)])
npts = 1
for c in choices:
    npts *= len(c)
print("number of monic points on N:", npts)
pts = [prod(c) for c in itertools.product(*choices)]
for z in pts[:5]:
    assert mul(z, conj(z)) == N4
fam = set()
for (w, k, c) in forms_and_sides(d):
    for S in (k, c):
        if 0 in S:
            z = prod(linfac(T, v) for v in sorted(S) if v != 0)
        else:
            z = prod(linfac(T, v) for v in sorted(S) if v != 1)
        fam.add(z)
# family points on this parent: F_S / (Y) or F_S/(Y+1)
famidx = [i for i, z in enumerate(pts) if z in fam]
print("Astra family points found on N:", len(famidx))
G = nx.Graph()
G.add_nodes_from(range(len(pts)))
for i, j in itertools.combinations(range(len(pts)), 2):
    if 2 * deg(add(pts[i], pts[j])) <= D:
        G.add_edge(i, j)
cl, w = nx.max_weight_clique(G, weight=None)
print("max c=0 cluster on N (all points):", len(cl), "; contains family points:", len(set(cl) & set(famidx)))
cl2, _ = nx.max_weight_clique(G.subgraph(famidx), weight=None)
print("max c=0 cluster among family points on N:", len(cl2))
