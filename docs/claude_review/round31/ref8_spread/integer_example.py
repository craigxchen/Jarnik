"""Referee's independent recomputation of the ten-mark example N = 1176852625 (Part VIII, COMPUTED).
Enumerates all representations of N, finds the cluster, normalises, and computes affine ranks
(exponent rows; plus residue indicators mod 3; plus mod 7), Smith invariants, coordinate orders,
windings and the exact phase.  Exact arithmetic (sympy Rationals, integers); mpmath only to round
the windings, with a margin check."""
from math import isqrt, gcd
from itertools import combinations
import mpmath as mp
from sympy import Matrix, Rational, factorint, ZZ, lcm
from sympy.matrices.normalforms import smith_normal_form

mp.mp.dps = 60
N = 1176852625
print("factorisation", factorint(N))
reps = []
for x in range(0, isqrt(N) + 1):
    y2 = N - x * x
    y = isqrt(y2)
    if y * y == y2:
        reps.append((x, y))
# all points in the first quadrant (x>0, y>=0)
pts = sorted({(x, y) for x, y in reps if x > 0}, key=lambda p: mp.atan2(p[1], p[0]))
print("points in first quadrant:", len(pts))
R = mp.sqrt(N)
sqR = mp.sqrt(R)
# best 10-point window by arc length over the full circle (use symmetry: angles mod 2pi)
allpts = []
for (x, y) in pts:
    for (a, b) in ((x, y), (-y, x), (-x, -y), (y, -x)):
        allpts.append((a, b))
allpts = sorted(set(allpts), key=lambda p: mp.atan2(p[1], p[0]) % (2 * mp.pi))
ang = [mp.atan2(p[1], p[0]) % (2 * mp.pi) for p in allpts]
n = len(allpts)
best = None
for i in range(n):
    j = (i + 9) % n
    da = (ang[j] - ang[i]) % (2 * mp.pi)
    L = R * da
    if best is None or L < best[0]:
        best = (L, i)
L, i0 = best
marks = [allpts[(i0 + k) % n] for k in range(10)]
print("best 10-point arc: C =", mp.nstr(L / sqR, 12), "marks:", marks)
assert L / sqR < mp.mpf(33) / 4

# Gaussian arithmetic
def mul(z, w): return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])
def conj(z): return (z[0], -z[1])
def nrm(z): return z[0]**2 + z[1]**2
def divides(d, z):
    a, b = mul(z, conj(d)); m = nrm(d)
    return a % m == 0 and b % m == 0
def gdiv(z, d):
    a, b = mul(z, conj(d)); m = nrm(d)
    return (a // m, b // m)
def ggcd(z, w):
    while w != (0, 0):
        # rounded division
        a, b = mul(z, conj(w)); m = nrm(w)
        qq = ((2*a + m) // (2*m), (2*b + m) // (2*m))
        r = (z[0] - mul(qq, w)[0], z[1] - mul(qq, w)[1])
        z, w = w, r
    return z
g = marks[0]
for z in marks[1:]:
    g = ggcd(g, z)
print("common Gaussian gcd norm:", nrm(g))
assert nrm(g) == 1
primes = {}
for p in factorint(N):
    assert p % 4 == 1
    for a in range(1, isqrt(p) + 1):
        b = isqrt(p - a * a)
        if a * a + b * b == p and a > b:
            primes[p] = (a, b); break
print("chosen pi_p:", primes)
plist = sorted(primes)
e = [factorint(N)[p] for p in plist]
def val(z, pi):
    c = 0
    while divides(pi, z):
        z = gdiv(z, pi); c += 1
    return c
A = [[val(z, primes[p]) for p in plist] for z in marks]
print("exponent rows a_ip:", A, "e =", e)
for col, ep in zip(zip(*A), e):
    assert min(col) == 0 and max(col) == ep
def affrank(rows):
    return Matrix([[r[k] - rows[0][k] for k in range(len(rows[0]))] for r in rows[1:]]).rank()
print("affine rank of exponent rows:", affrank(A))
feat = [list(r) for r in A]
for ell in (3, 7, 11, 19):
    cls = [(x % ell, y % ell) for x, y in marks]
    for z in cls:
        assert (z[0]**2 + z[1]**2 - N) % ell == 0
    used = sorted(set(cls))
    for r, c in zip(feat, cls):
        r.extend(int(c == u) for u in used)
    print(f"  after residue indicators mod {ell} ({len(used)} classes used of {ell+1}): affine rank {affrank(feat)}")

B = Matrix([[2 * (A[i][k] - A[0][k]) for k in range(5)] for i in range(1, 10)])
snf = smith_normal_form(B, domain=ZZ)
print("Smith invariants of B:", [abs(snf[k, k]) for k in range(5)])
# orders of e_p modulo the row lattice L (full rank 5): smallest D with D e_p in L.
# Solve x B = D e_p : x = D e_p (B^T B)^{-1} B^T restricted ... use rational solve of B^T x = e_p
orders = []
for k in range(5):
    ek = Matrix([1 if t == k else 0 for t in range(5)])
    # e_p in L_Q: find rational y with B^T y = e_p; minimal D s.t. some integral y exists: use HNF of rows
    from sympy.matrices.normalforms import hermite_normal_form
    H = hermite_normal_form(B.T)  # columns form basis of column lattice of B^T = row lattice of B
    # basis vectors = columns of H (5x5)
    coeff = H.inv() * ek
    orders.append(int(lcm([Rational(c).q for c in coeff])))
print("orders of coordinate vectors mod L:", orders)
# windings
t = [mp.atan2(primes[p][1], primes[p][0]) / (mp.pi / 2) for p in plist]
Bt = [sum(int(B[i, k]) * t[k] for k in range(5)) for i in range(9)]
kw = [int(mp.nint(x)) for x in Bt]
us = [x - k for x, k in zip(Bt, kw)]
print("windings k:", kw, " max|u| =", mp.nstr(max(abs(x) for x in us), 5))
assert max(abs(x) for x in us) < 0.25  # arc width in quarter turns is 0.0283
sol = B.gauss_jordan_solve(Matrix(kw))
tstar = sol[0]
print("exact phase t* (quarter turns):", list(tstar), " free params:", sol[1].shape)
assert B * tstar == Matrix(kw)
print("t - t* :", [mp.nstr(t[k] - float(tstar[k]), 5) for k in range(5)])
