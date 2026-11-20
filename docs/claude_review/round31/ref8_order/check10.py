# Referee check of VIII.2 on the 10-point record cluster N = 1176852625.
# Exact integer / rational arithmetic; mpmath only for the angles.
from math import isqrt, gcd
from fractions import Fraction as F
import sympy
from sympy import Matrix, ZZ, lcm
from sympy.matrices.normalforms import smith_normal_form, hermite_normal_form
import mpmath as mp
mp.mp.dps = 60

N = 1176852625
fac = sympy.factorint(N)
print("N =", N, fac)

# all lattice points
pts = []
for x in range(-isqrt(N), isqrt(N) + 1):
    y2 = N - x * x
    y = isqrt(y2)
    if y * y == y2:
        pts.append((x, y))
        if y: pts.append((x, -y))
print("number of lattice points", len(pts))
R = mp.sqrt(N)
ang = sorted(pts, key=lambda p: mp.atan2(p[1], p[0]))
th = [mp.atan2(p[1], p[0]) for p in ang]
n = len(ang)
best = None
for i in range(n):
    j = i + 9
    d = th[j % n] - th[i] + (2 * mp.pi if j >= n else 0)
    C = d * R / mp.sqrt(R)
    if best is None or C < best[0]:
        best = (C, i)
C, i0 = best
cl = [ang[(i0 + k) % n] for k in range(10)]
print("best 10-window arc constant C =", mp.nstr(C, 12))
print("cluster", cl)

# Gaussian arithmetic
def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def conj(a): return (a[0], -a[1])
def gnorm(a): return a[0]*a[0]+a[1]*a[1]
def gdivexact(a, b):
    nb = gnorm(b); num = gmul(a, conj(b))
    if num[0] % nb or num[1] % nb: return None
    return (num[0]//nb, num[1]//nb)
def gmod_round(a, b):
    nb = gnorm(b); num = gmul(a, conj(b))
    q = ((2*num[0]+nb)//(2*nb), (2*num[1]+nb)//(2*nb))
    qb = gmul(q, b); return (a[0]-qb[0], a[1]-qb[1])
def ggcd(a, b):
    while b != (0, 0):
        a, b = b, gmod_round(a, b)
    return a
def val(z, pi):
    k = 0
    while True:
        q = gdivexact(z, pi)
        if q is None: return k
        z = q; k += 1

g = cl[0]
for z in cl[1:]: g = ggcd(g, z)
print("common gcd", g, "norm", gnorm(g), sympy.factorint(gnorm(g)))
zs = [gdivexact(z, g) for z in cl]
Np = gnorm(zs[0]); assert all(gnorm(z) == Np for z in zs)
print("N' =", Np, sympy.factorint(Np))
W = mp.log(Np)
Rp = mp.sqrt(Np)
# arc constant after normalisation
thz = [mp.atan2(z[1], z[0]) for z in zs]
# primes: choose pi_p = a+bi with a>b>0
primes = sorted(sympy.factorint(Np))
pis = {}
for p in primes:
    for a in range(1, isqrt(p)+1):
        b2 = p - a*a; b = isqrt(b2)
        if b*b == b2 and a > b > 0: pis[p] = (a, b)
ep = sympy.factorint(Np)
print("primes", primes, "pi_p", pis, "e_p", ep)
A = [[val(z, pis[p]) for p in primes] for z in zs]
for row, z in zip(A, zs):
    # check v_pi + v_pibar = e_p
    assert all(val(z, conj(pis[p])) == ep[p] - row[k] for k, p in enumerate(primes))
print("exponent rows a_i:", A)
M = len(zs); r = len(primes)
Bm = Matrix([[2*(A[i][k]-A[0][k]) for k in range(r)] for i in range(1, M)])
s = Bm.rank()
print("rank B =", s, " r =", r)
snf = smith_normal_form(Bm, domain=ZZ)
print("Smith invariants of B:", [snf[k, k] for k in range(min(snf.shape)) if snf[k, k] != 0])
# row lattice basis via HNF of B^T (column-style) -> take nonzero rows of HNF of B
H = hermite_normal_form(Bm.T).T  # rows generate same lattice? verify
# safer: basis from LLL-free approach: use sympy's rref over Z through Smith: compute lattice basis by HNF of Bm (row HNF)
def row_hnf_basis(Mx):
    Hc = hermite_normal_form(Mx.T)  # columns generate column lattice of Mx^T = row lattice of Mx
    return Hc.T
Bas = row_hnf_basis(Bm)
Bas = Matrix([list(Bas.row(k)) for k in range(Bas.rows) if any(Bas.row(k))])
print("row-lattice basis (rows):", Bas.tolist(), "det", Bas.det())
# check Bas rows lie in L and B rows lie in span_Z(Bas)
Binv = Bas.inv()
for k in range(Bm.rows):
    x = Bm.row(k) * Binv
    assert all(xx.q == 1 for xx in x)
for k in range(Bas.rows):
    # Bas row in L: solve y B = Bas_k integrally? check via determinant equality of lattices
    pass
detL = abs(Bas.det())
print("index [Z^r : L] =", detL, " product of Smith invariants =", sympy.prod([snf[k, k] for k in range(min(snf.shape)) if snf[k, k] != 0]))
orders = []
for k in range(r):
    e = Matrix([[1 if j == k else 0 for j in range(r)]])
    x = e * Binv
    orders.append(int(lcm([xx.q for xx in x])))
print("coordinate orders D_L(e_p):", dict(zip(primes, orders)))

# phases: t_p = 2 arg(pi_p)/pi  (quarter-turns)
t = [2*mp.atan2(pis[p][1], pis[p][0])/mp.pi for p in primes]
# u_i = (2/pi)(theta_i - theta_1), lifted to (-pi, pi]
def wrap(x):
    while x > mp.pi: x -= 2*mp.pi
    while x <= -mp.pi: x += 2*mp.pi
    return x
u_true = [2*wrap(thz[i]-thz[0])/mp.pi for i in range(1, M)]
Bt = [sum(Bm[i, k]*t[k] for k in range(r)) for i in range(M-1)]
kvec = [int(mp.nint(Bt[i]-u_true[i])) for i in range(M-1)]
u = [Bt[i]-kvec[i] for i in range(M-1)]
# u should equal u_true modulo integers -> check equality mod 1 handled by kvec choice
print("max |u_i - (2/pi)(theta_i-theta_1)| =", mp.nstr(max(abs(u[i]-u_true[i]) for i in range(M-1)), 5))
uinf = max(abs(x) for x in u)
Cn = (thz[-1]-thz[0]) if True else 0
Cnorm = abs(wrap(max(thz)-min(thz))) * mp.sqrt(Rp)
print("normalised arc constant C' =", mp.nstr(Cnorm, 10), " ||u||_inf =", mp.nstr(uinf, 8),
      " bound (2C'/pi) e^{-W/4} =", mp.nstr(2*Cnorm/mp.pi*mp.e**(-W/4), 8))
# max minor H
from itertools import combinations
Hmax = 1
for j in range(1, s+1):
    for I in combinations(range(M-1), j):
        for J in combinations(range(r), j):
            d = abs(Bm.extract(list(I), list(J)).det())
            if d > Hmax: Hmax = d
print("H (max |minor|) =", Hmax, " log H / W =", mp.nstr(mp.log(Hmax)/W, 6), " W =", mp.nstr(W, 8))
print("(s+1) H ||u|| =", mp.nstr((s+1)*Hmax*uinf, 8), " (needs < 1)")
# choose nonsingular s-minor: pick the one maximising |det| (any works if condition holds)
bestI = None
for I in combinations(range(M-1), s):
    for J in combinations(range(r), s):
        d = abs(Bm.extract(list(I), list(J)).det())
        if d and (bestI is None or d < bestI[0]):
            bestI = (d, I, J)
d, I, J = bestI
T = Bm.extract(list(I), list(J))
print("chosen minor rows", I, "cols", J, "d =", d)
# exact t*: B_I t* = k_I on cols J, other cols keep t (here J = all cols when s = r)
# when s == r, t* = T^{-1} k_I exactly rational
assert s == r
tstar = T.inv() * Matrix([kvec[i] for i in I])
print("t* (quarter-turns, exact):", [str(x) for x in tstar])
print("t* mod 1:", [str(F(int(x.p), int(x.q)) % 1) for x in tstar])
ok = all(sum(Bm[i, k]*tstar[k] for k in range(r)) == kvec[i] for i in range(M-1))
print("B t* = k on ALL rows:", ok)
print("|t - t*|_inf =", mp.nstr(max(abs(t[k]-mp.mpf(tstar[k].p)/tstar[k].q) for k in range(r)), 6))
print("denominators of t*_p:", [tstar[k].q for k in range(r)], " vs orders", orders)
# order-height numbers
for k, p in enumerate(primes):
    D = orders[k]
    dist = abs(D*t[k] - mp.nint(D*t[k]))
    alpha_pow = abs(mp.e**(2j*mp.pi*D*t[k]) - 1)
    print(f"p={p}: D={D}, ||D t_p|| = {mp.nstr(dist,6)}, -log||D t_p|| = {mp.nstr(-mp.log(dist),5)}, "
          f"(1+log p) log(2D) = {mp.nstr((1+mp.log(p))*mp.log(2*D),5)}, ratio to W = {mp.nstr((1+mp.log(p))*mp.log(2*D)/W,4)}")
