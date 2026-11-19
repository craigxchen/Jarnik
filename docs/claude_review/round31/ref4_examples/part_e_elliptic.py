"""(e) A 12-mark configuration over Q with delta = 0 and all E_ij = 0.

E: y^2 = x(x+r^2)(x+s^2) has torsion containing G = Z/2 x Z/4.  Embed E in
P^1 x P^1 as a (2,2)-curve by two degree-2 maps f_a, f_b (pencils of lines
through points of E).  T = P0 + G (P0 of infinite order) gives 8 points of E in
8 distinct t-fibres.  A (1,1)-curve (graph of a Moebius map: column
(-(b t + d), a t + c)) cuts E in 4 points; any 4 points of E with the right sum
lie on one.  Every 4-subset of T with the right sum gives a mark: alpha_i = 4,
L = 8, deg D_ij <= 2 = (4+4)/2 - 8/4.  Layers: the 8 points of T."""
import itertools, sys
from fractions import Fraction as Fr
import sympy as sp
from front import t, contact_data, iv2_data, key
from engine import chern_side
from collections import Counter

r, s = int(sys.argv[1]), int(sys.argv[2])
B, C = r * r + s * s, r * r * s * s      # y^2 = x^3 + B x^2 + C x
O = None


def on(P):
    return P is None or P[1] ** 2 == P[0] ** 3 + B * P[0] ** 2 + C * P[0]


def add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and y1 == -y2:
        return None
    if P == Q:
        lam = (3 * x1 * x1 + 2 * B * x1 + C) / (2 * y1)
    else:
        lam = (y2 - y1) / (x2 - x1)
    x3 = lam * lam - B - x1 - x2
    y3 = lam * (x1 - x3) - y1
    return (x3, y3)


def neg(P):
    return None if P is None else (P[0], -P[1])


def mulp(n, P):
    R = None
    Q = P if n >= 0 else neg(P)
    for _ in range(abs(n)):
        R = add(R, Q)
    return R


def find_point():
    from math import isqrt
    for v in range(1, 10):
        for u in range(-80 * v * v, 80 * v * v + 1):
            x = Fr(u, v * v)
            rhs = x * (x * x + B * x + C)
            if rhs <= 0:
                continue
            n, d = rhs.numerator, rhs.denominator
            if isqrt(n) ** 2 == n and isqrt(d) ** 2 == d:
                P = (x, Fr(isqrt(n), isqrt(d)))
                # non-torsion: Mazur, order <= 16 checked
                Q = P
                ok = True
                for k in range(2, 17):
                    Q = add(Q, P)
                    if Q is None:
                        ok = False
                        break
                if ok:
                    return P
    raise SystemExit("no point")


P0 = find_point()
r_, s_ = Fr(r), Fr(s)
# order-4 point: x = r s, y^2 = rs (rs + r^2)(rs + s^2) = r^2 s^2 (r+s)^2
P4 = (r_ * s_, r_ * s_ * (r_ + s_))
assert on(P4) and mulp(4, P4) is None and mulp(2, P4) == (Fr(0), Fr(0))
T2 = (Fr(-r * r), Fr(0))
G = [add(mulp(i, P4), mulp(j, T2)) for i in range(4) for j in range(2)]
assert len(set(G)) == 8 and all(on(g) for g in G)
T = [add(P0, g) for g in G]
target = add(mulp(4, P0), (Fr(0), Fr(0)))      # the best sum class (0,2)
sig_a = neg(P0)
sig_b = add(target, neg(sig_a))


def slope_from(Abase, P):
    # slope of the line through Abase and P (tangent if equal)
    if P == Abase:
        return (3 * P[0] ** 2 + 2 * B * P[0] + C) / (2 * P[1])
    assert P[0] != Abase[0], "pole"
    return (P[1] - Abase[1]) / (P[0] - Abase[0])


Aa, Ab = neg(sig_a), neg(sig_b)
pts = [(slope_from(Aa, P), slope_from(Ab, P)) for P in T]
assert len(set(p[0] for p in pts)) == 8, "t-fibres not distinct"


def det4(M):
    return sp.Matrix(M).det()


def mobius(S):
    rows = [[x * y, x, y, 1] for x, y in (pts[k] for k in S)]
    M = sp.Matrix(rows)
    if M.det() != 0:
        return None
    ns = M[:3, :].nullspace()
    assert len(ns) == 1
    a, b, c, d = [sp.nsimplify(z) for z in ns[0]]
    return a, b, c, d


sets, cols = [], []
allcnt = 0
for S in itertools.combinations(range(8), 4):
    mb = mobius(S)
    if mb is None:
        continue
    allcnt += 1
    a, b, c, d = mb
    assert b * c - a * d != 0
    sets.append(S)
    cols.append((-(b * t + d), a * t + c))
sumcls = [S for S in itertools.combinations(range(8), 4)
          if add(add(T[S[0]], T[S[1]]), add(T[S[2]], T[S[3]])) == target]
print("P0 =", P0, "; 4-subsets on (1,1)-curves:", allcnt, "; 4-subsets in the sum class:", len(sumcls))
assert sorted(sumcls) == sorted(sets)
m = len(sets)
cd = contact_data(cols)
N, a = {}, [dict() for _ in range(m)]
for k, (x, y) in enumerate(pts):
    p = sp.Poly(t - sp.Rational(x.numerator, x.denominator), t, domain='QQ')
    kk = key(p)
    N[kk] = (p, 1)
    for i, S in enumerate(sets):
        if k in S:
            a[i][kk] = 1
ar = iv2_data(cd, N, a)
ch = chern_side(m, cd['finite'], cd['vinf'], cd['degs'])
print("s_h =", [sum(1 for S in sets if k in S) for k in range(8)])
print("deg D_ij multiset:", Counter(d.degree() for d in cd['D'].values()))
print(f"m={m} L={ar['L']} delta={ar['delta']} sumE={ar['sumE']} IV.2 lhs={ar['lhs']} rhs={ar['rhs']} slack={ar['slack']}")
print(f"Chern side: S={ch['S']} I={ch['I']} big={ch['big']} chains={ch['chains']} unblown={ch['unblown']} "
      f"3e_orb-(K+B)^2={ch['Q']} ledger={ch['ledger']} bound={ch['bound']}")
