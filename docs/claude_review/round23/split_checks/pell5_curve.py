"""Exact curve data of the m-czech Pell family as a curve in split LM_5.
Function field K = Q(x)[z]/(z^2-56x^2-65) (a conic with rational point (1,11), so K = Q(t)).
Finite places of K over x=c split (56c^2+65 square) -> place (c, z0) with z0 = +-sqrt(56c^2+65).
Valuation of a divisor f(x,z) at place (c,z0): ord_{x=c} of f(x, z) restricted... we compute via norm: f is a product of
places; test f(c, z0) == 0. Since all 16 places are simple zeros of n (x^2-c^2 simple), ord is 0/1.
"""
import itertools, sympy as sp
x, z = sp.symbols('x z')
cs = [1, -1, 2, -2, 4, -4, 7, -7]
places = []
for c in cs:
    r = sp.sqrt(56*c*c+65); assert r.is_integer
    places += [(c, r), (c, -r)]
# candidate divisors in window, identified numerically at a large Pell solution
Xn, Zn = 803080430161, None
import math
Zn = math.isqrt(56*Xn*Xn+65); assert Zn*Zn == 56*Xn*Xn+65
n = (Xn**2-1)*(Xn**2-4)*(Xn**2-16)*(Xn**2-49)
s = math.isqrt(n); q4 = n**0.25
lin = {c: x + c for c in cs}  # factor (x+c) vanishes at x=-c
cands = {}
for r in range(9):
    for comb in itertools.combinations(cs, r):
        f = sp.Integer(1)
        for c in comb: f *= (x + c)
        cands[('lin',) + comb] = sp.expand(f)
cands[('f5',)] = x**4 - 35*x**2 + 56 + 2*x*z
win = []
for key, f in cands.items():
    v = int(f.subs({x: Xn, z: Zn}))
    if v > 0 and n % v == 0 and s < v <= s + 40*q4 + 1:
        win.append((v, key, f))
win.sort()
print("divisors in window:")
for v, key, f in win: print("  offset/n^(1/4) = %.4f" % ((v - s)/q4), key, sp.factor(f) if key[0]=='lin' else f)
F = [f for v, key, f in win]
# valuation vectors: f vanishes at place (c,z0)?  for x+c' factor: vanishes at x=-c' both z
V = []
for (c, z0) in places:
    vec = tuple(1 if sp.simplify(f.subs({x: c, z: z0})) == 0 else 0 for f in F)
    V.append(((c, z0), vec))
BC = 0
cutcount = {}
for p, vec in V:
    w = max(vec) - min(vec); BC += w
    A = frozenset(j for j in range(5) if vec[j] == 1)
    key = min(A, frozenset(range(5)) - A, key=lambda S: sorted(S))
    if w: cutcount[tuple(sorted(key))] = cutcount.get(tuple(sorted(key)), 0) + 1
    print("  place", p, "vec", vec)
print("B.C (finite places; infinite places contribute 0 since all degrees 4) =", BC)
print("cuts hit:", len(cutcount), cutcount)
# order at infinity: d_j - d_1 as function; leading degree in x of differences (z ~ sqrt56 x)
t = sp.symbols('t')
for j in range(1, 5):
    diff = sp.expand(F[j] - F[0])
    print("  F%d-F0 =" % j, diff)
