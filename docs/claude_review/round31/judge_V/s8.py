# lower's by-product: an 8-point c=0 cluster inside the family at d=4 with point degree 104 (< 112 of V.1).
import itertools
from f4 import *
# F_4 elements 0,1,2=theta,3=theta^2 ; mult table
ML = [[0,0,0,0],[0,1,2,3],[0,2,3,1],[0,3,1,2]]
INV = [None, 1, 3, 2]
def tolist(x):
    A, B = x; n = max(A.bit_length(), B.bit_length())
    return [((A >> i) & 1) | (((B >> i) & 1) << 1) for i in range(n)]
def trim(p):
    while p and p[-1] == 0: p.pop()
    return p
def pmod(a, b):
    a = a[:]; db = len(b) - 1; ib = INV[b[-1]]
    while len(a) - 1 >= db and a:
        c = ML[a[-1]][ib]; s = len(a) - 1 - db
        for i, bi in enumerate(b): a[s + i] ^= ML[c][bi]
        trim(a)
    return a
def pgcd(a, b):
    while b: a, b = b, pmod(a, b)
    return a
def par(x): return bin(x).count('1') & 1
d = 4; Q = 16; L = 15; forms = [c for c in range(Q) if c & 1]
lin = lambda v: ((1 << L) ^ v, 1)
best = None
for sides in itertools.product((0, 1), repeat=8):
    F = []
    for c, s in zip(forms, sides):
        p = ONE
        for v in range(Q):
            if par(c & v) == s: p = mul4(p, lin(v))
        F.append(p)
    D0 = deg4(F[0])
    G = tolist(F[0])
    for f in F[1:]: G = pgcd(G, tolist(f))
    g = len(G) - 1
    md = max(deg4(add4(a, b)) for a, b in itertools.combinations(F, 2))
    if 2 * (md - g) <= D0 - g:
        if best is None or D0 - g < best[0]: best = (D0 - g, sides, g, md - g)
print("d=4, T=t^15, all 8 forms, best reduced point degree, sides, gcd deg, max reduced diff:", best)
