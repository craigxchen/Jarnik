"""(e) Kummer 16_6 configuration on a split quadric = 16 Moebius sections with
delta = 0 and all E_ij = 0 at L = 16 (each section through 6 of the 16 layer
points, any two sharing exactly 2 of them: a biplane of order 4).

Heisenberg orbit of p in P^3 (coordinates indexed by F_2^2: 00,01,10,11): the
16 points (signs) x (double transpositions) of p.  If p0^2+p1^2 = p2^2+p3^2 the
whole orbit lies on Q: x0^2+x1^2 = x2^2+x3^2, i.e. X Y' = Z W with
X = x0-x2, Y' = x0+x2, Z = x3-x1, W = x3+x1; P^1 x P^1 coordinates
t = W/X, y = Z/X.  Plane sections a X + b Y' + c Z + d W = 0 are the graphs
y = -(a + d t)/(b t + c)."""
import itertools, sys
import sympy as sp
from fractions import Fraction as Fr
from front import t, contact_data, iv2_data, key
from engine import chern_side
from collections import Counter

p = [int(z) for z in sys.argv[1:5]] if len(sys.argv) >= 5 else [1, 8, 4, 7]
assert p[0] ** 2 + p[1] ** 2 == p[2] ** 2 + p[3] ** 2
idx = [(0, 0), (0, 1), (1, 0), (1, 1)]
pts3 = set()
for a in idx:
    for b in idx:
        q = [0] * 4
        for k, v in enumerate(idx):
            w = ((v[0] + a[0]) % 2, (v[1] + a[1]) % 2)
            sgn = (-1) ** ((b[0] * v[0] + b[1] * v[1]) % 2)
            q[k] = sgn * p[idx.index(w)]
        # projective normalisation
        g = next(z for z in q if z != 0)
        pts3.add(tuple(Fr(z, g) for z in q))
pts3 = sorted(pts3)
assert len(pts3) == 16
for x in pts3:
    assert x[0] ** 2 + x[1] ** 2 == x[2] ** 2 + x[3] ** 2
# planes with >= 6 points
M = sp.Matrix
planes = {}
for i, j, k in itertools.combinations(range(16), 3):
    A = M([list(pts3[i]), list(pts3[j]), list(pts3[k])])
    ns = A.nullspace()
    if len(ns) != 1:
        continue
    pl = ns[0]
    pl = pl / next(z for z in pl if z != 0)
    on = tuple(l for l in range(16) if (M([list(pts3[l])]) * pl)[0] == 0)
    if len(on) >= 6:
        planes[tuple(pl)] = on
print("planes through >= 6 of the 16 points:", len(planes), Counter(len(v) for v in planes.values()))
blocks = list(planes.values())
print("block intersections:", Counter(len(set(a) & set(b)) for a, b in itertools.combinations(blocks, 2)))


def tcoord(x):
    X, Yp, Z, W = x[0] - x[2], x[0] + x[2], x[3] - x[1], x[3] + x[1]
    assert X * Yp == Z * W
    if X != 0:
        return (W / X, Z / X)
    raise ValueError("point on X=0; choose another p")


TY = [tcoord(x) for x in pts3]
print("distinct t-coordinates:", len(set(z[0] for z in TY)), " distinct y:", len(set(z[1] for z in TY)))
assert len(set(z[0] for z in TY)) == 16
cols = []
for pl, on in planes.items():
    c0, c1, c2, c3 = [sp.Rational(z.p, z.q) if isinstance(z, sp.Rational) else sp.nsimplify(z) for z in pl]
    # plane c . x = 0 with x0 = (X+Y')/2, x2 = (Y'-X)/2, x1 = (W-Z)/2, x3 = (W+Z)/2
    a_ = (c0 - c2) / 2      # coefficient of X
    b_ = (c0 + c2) / 2      # of Y'
    cZ = (c3 - c1) / 2      # of Z
    dW = (c1 + c3) / 2      # of W
    cols.append((sp.expand(-(a_ + dW * t)), sp.expand(b_ * t + cZ)))
    for l in on:
        tt, yy = TY[l]
        assert a_ + b_ * tt * yy + cZ * yy + dW * tt == 0
m = len(cols)
cd = contact_data(cols)
N, a = {}, [dict() for _ in range(m)]
for l, (tt, yy) in enumerate(TY):
    pp = sp.Poly(t - sp.Rational(tt.numerator, tt.denominator), t, domain='QQ')
    N[key(pp)] = (pp, 1)
    for i, on in enumerate(blocks):
        if l in on:
            a[i][key(pp)] = 1
ar = iv2_data(cd, N, a)
ch = chern_side(m, cd['finite'], cd['vinf'], cd['degs'])
print("deg D_ij:", Counter(d.degree() for d in cd['D'].values()), " s_h:", Counter(sum(1 for on in blocks if l in on) for l in range(16)))
print(f"m={m} L={ar['L']} delta={ar['delta']} sumE={ar['sumE']} IV.2 lhs={ar['lhs']} rhs={ar['rhs']} slack={ar['slack']}")
print(f"Chern: S={ch['S']} I={ch['I']} infsteps={ch['steps']} big={Counter(ch['big'])} chains={ch['chains']} unblown={ch['unblown']} "
      f"(K+B)^2={ch['KB2']} e_orb={ch['eorb']} 3e-(K+B)^2={ch['Q']} ledger={ch['ledger']} bound={ch['bound']}")

# extra sections of class (1,d), d = 0..3, through >= 2d + L/4 = 2d + 4 of the points
import numpy as np
def extra(d, need):
    found = []
    mons = [(i, e) for e in (0, 1) for i in range(d + 1)]   # y^e t^i
    for S in itertools.combinations(range(16), need):
        A = M([[TY[l][0] ** i * TY[l][1] ** e for i, e in mons] for l in S])
        if A.rank() < len(mons):
            found.append(S)
    return found
if '--extra' in sys.argv:
    for d, need in ((0, 4), (1, 7), (2, 8)):
        f = extra(d, need)
        print(f"curves of class (1,{d}) through >= {need} of the 16 points: {len(f)} subsets")
