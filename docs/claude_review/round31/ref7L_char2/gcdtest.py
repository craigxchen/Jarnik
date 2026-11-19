# (f) sanity test of the VII.2 gcd estimate in characteristic 2: fixed-box coprime F1,F2 over F_2,
# common vanishing deg gcd(F1(t^n),F2(t^n)) for n=(L, aL+b), L up to 2048 (incl. powers of 2).
import random, itertools
from gf import *
random.seed(5)
def spec(F, n):
    deg_ = max(e[0]*n[0] + e[1]*n[1] for e in F); r = [0]*(deg_+1)
    for e in F: r[e[0]*n[0]+e[1]*n[1]] ^= 1
    return p_trim(r)
def strip(p):  # remove powers of t
    i = 0
    while i < len(p) and p[i] == 0: i += 1
    return p[i:]
box = [(i, j) for i in range(3) for j in range(3)]
pairs = []
while len(pairs) < 6:
    F1 = set(random.sample(box, 4)); F2 = set(random.sample(box, 4))
    # coprimality in F_2[Z1,Z2] tested crudely via specialisation at a generic-ish n with coprime result
    if F1 != F2: pairs.append((F1, F2))
# add a char-2 flavoured pair: 1+Z1+Z2 and 1+Z1^2 Z2 + Z2^2
pairs.append(({(0,0),(1,0),(0,1)}, {(0,0),(2,1),(0,2)}))
for F1, F2 in pairs:
    row = []
    for L in (64, 128, 256, 512, 1024, 2048):
        best = 0
        for (a, b) in ((1, 1), (2, 1), (3, 2), (1, 3)):
            n = (L, a*L + b)
            g = p_gcd(strip(spec(F1, n)), strip(spec(F2, n)))
            best = max(best, deg(g))
        row.append(best)
    print(sorted(F1), sorted(F2), "max gcd degree over directions, L=64..2048:", row)
