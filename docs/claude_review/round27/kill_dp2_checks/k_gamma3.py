# Killer: (i) sections of -dK with ord_P >= a, ord_P' >= b (dimension), Gamma_k classes for k = 2,3,4 and -2;
# (ii) D-points from Dirichlet approximants of the real cubic slope; check x + nu_+ against Gamma_3 bound 3/4.
import sympy as sp, math
from fractions import Fraction as Fr
x, y, z, e, s, W = sp.symbols('x y z e s W')
f = x**3*y + 2*x**2*y*z - 2*x*y**3 - 2*x*y*z**2 + 2*x*z**3 + 2*y**4 + 2*y**3*z + y**2*z**2 - y*z**3 + z**4
g = sp.expand(f.subs({z: 1}).subs({x: e*x, y: e*y}))
ORD = 16
wser = sp.expand(sp.series(sp.sqrt(g), e, 0, ORD).removeO())

def dim_sections(d, a, b):
    mons = [x**i*y**j*z**(d-i-j) for i in range(d+1) for j in range(d+1-i)]
    mons += [s*x**i*y**j*z**(d-2-i-j) for i in range(d-1) for j in range(d-1-i)]
    rows = []
    for sign, order in ((1, a), (-1, b)):
        if order <= 0: continue
        cols = []
        for mo in mons:
            lo = sp.expand(mo.subs({s: W}).subs({z: 1}).subs({x: e*x, y: e*y}).subs({W: sign*wser}))
            col = {}
            for dd in range(order):
                cd = sp.expand(lo.coeff(e, dd))
                if cd != 0:
                    for mon, c in sp.Poly(cd, x, y).terms():
                        col[(dd, mon)] = c
            cols.append(col)
        keys = sorted(set().union(*cols))
        rows += [[col.get(k, 0) for col in cols] for k in keys]
    M = sp.Matrix(rows) if rows else sp.zeros(0, len(mons))
    return len(mons), len(mons) - (M.rank() if rows else 0)

for (d, a, b) in [(2, 3, 0), (2, 0, 3), (6, 8, 3), (6, 9, 3), (6, 8, 4), (6, 3, 8), (4, 5, 1), (4, 6, 0), (3, 4, 1), (3, 5, 0)]:
    n, dim = dim_sections(d, a, b)
    print(f"d={d} ordP>={a} ordP'>={b}: h0(-dK)={n}, dim={dim}", flush=True)

# D-points: line of slope (a',b') meets pi(D) again at (-8G3 a' : -8G3 b' : G4)
a1, b1 = sp.symbols('a1 b1')
G3 = sp.Poly(8*(a1**3 + sp.Rational(5, 2)*a1**2*b1 - sp.Rational(5, 4)*a1*b1**2 + sp.Rational(19, 8)*b1**3), a1, b1)
G4 = sp.Poly(64*(-a1**4/4 + a1**3*b1/2 + a1**2*b1**2/8 - sp.Rational(13, 8)*a1*b1**3 + sp.Rational(119, 64)*b1**4), a1, b1)
Tq = lambda A, B, C: Fr(-1, 2)*A*A - Fr(1, 2)*A*B + A*C + Fr(3, 8)*B*B - Fr(1, 2)*B*C + C*C
fz = lambda A, B, C: A**3*B + 2*A**2*B*C - 2*A*B**3 - 2*A*B*C**2 + 2*A*C**3 + 2*B**4 + 2*B**3*C + B**2*C**2 - B*C**3 + C**4
# real root theta of 19t^3-10t^2+20t+8 (t = b'/a'), continued fraction convergents
import mpmath
mpmath.mp.dps = 200
th = mpmath.findroot(lambda t: 19*t**3 - 10*t**2 + 20*t + 8, -0.33)
cf = []; v = th
for _ in range(40):
    q0 = int(mpmath.floor(v)); cf.append(q0); v = 1/(v - q0)
convs = []; h0, h1, k0, k1 = 1, cf[0], 0, 1
for q0 in cf[1:]:
    h0, h1 = h1, q0*h1 + h0; k0, k1 = k1, q0*k1 + k0
    convs.append((k1, h1))   # a' = denominator, b' = numerator (t = b'/a')
def oddpart(n):
    while n % 2 == 0: n //= 2
    return n
print("D-points (a'=k, b'=h):  log10 c, x, nu_+, x+nu_+  [Gamma_3 bound: x + nu_+ <= 3/4 + o(1)]")
for (ap, bp) in convs[2:16]:
    g3 = int(G3.eval({a1: ap, b1: bp})); g4 = int(G4.eval({a1: ap, b1: bp}))
    A, B, C = -8*g3*ap, -8*g3*bp, g4
    gg = math.gcd(math.gcd(A, B), C); A, B, C = A//gg, B//gg, C//gg
    if C < 0: A, B, C = -A, -B, -C
    sv = Tq(A, B, C)
    assert sv.denominator == 1 and sv*sv == fz(A, B, C)
    sv = int(sv)
    m = math.gcd(A, B)
    # every odd p | m has s == c^2 mod p (P-side) ?
    mo = oddpart(m)
    assert (sv - C*C) % mo == 0
    X = Fr(max(abs(A), abs(B)), C)
    lc = math.log(C)
    xx = -math.log(X)/lc; nup = math.log(mo)/lc
    print(f"  {lc/math.log(10):7.2f}  x={xx:.4f}  nu+={nup:.4f}  x+nu+={xx+nup:.4f}")
