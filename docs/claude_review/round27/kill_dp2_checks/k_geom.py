# Killer's independent geometry check for lane dp2 (own code).
import sympy as sp
from sympy import Rational as Q
x, y, z, t, a, b = sp.symbols('x y z t a b')
f = x**3*y + 2*x**2*y*z - 2*x*y**3 - 2*x*y*z**2 + 2*x*z**3 + 2*y**4 + 2*y**3*z + y**2*z**2 - y*z**3 + z**4
print("f(0,0,1) =", f.subs({x: 0, y: 0, z: 1}))

# 1. smoothness of plane quartic B: f=fx=fy=fz=0 has no point in P^2(Qbar)
fx, fy, fz = [sp.diff(f, v) for v in (x, y, z)]
def sing_check(mod=None):
    kw = dict(order='grevlex')
    if mod: kw['modulus'] = mod
    G = sp.groebner([fx, fy, fz], x, y, z, **kw)
    lms = [sp.Poly(g, x, y, z).monoms(order='grevlex')[0] for g in G.exprs]
    # pure powers of x, y, z among leading monomials => zero-dim ideal supported at origin => no projective point
    pp = {}
    for m in lms:
        nz = [i for i in range(3) if m[i] > 0]
        if len(nz) == 1:
            i = nz[0]; pp[i] = min(pp.get(i, 99), m[i])
    return pp
for mod in (None, 3, 5):
    print("pure-power leading monomials (var index: exponent), modulus", mod, ":", sing_check(mod))

# 2. local expansion at P: z=1, s = sqrt(f(x,y,1)), branch s(0,0)=1
e = sp.symbols('e')
g = sp.expand(f.subs({z: 1}).subs({x: e*x, y: e*y}))
s_ser = sp.series(sp.sqrt(g), e, 0, 6).removeO()
w = [sp.expand(s_ser.coeff(e, i)) for i in range(6)]
print("w1 =", w[1]); print("w2 =", w[2]); print("w3 =", w[3])
T = z**2 + z*w[1] + w[2]   # homogenised 2-jet (w1 linear, w2 quadratic in x,y)
print("T =", sp.expand(T))
rem = sp.expand(f - T**2)
F3 = sp.expand(rem.coeff(z, 1)); F4p = sp.expand(rem.coeff(z, 0))
print("f - T^2 has z-degree", sp.Poly(rem, z).degree(), "; F3 =", F3, "; F4' =", F4p)
print("2*w3 == F3 ?", sp.simplify(2*w[3] - F3) == 0)
# 3. resultant of binary forms F3, F4' (homogeneous), via dehomogenisation + leading checks
F3t = sp.Poly(F3.subs({x: 1, y: t}), t); F4t = sp.Poly(F4p.subs({x: 1, y: t}), t)
print("deg F3(1,t) =", F3t.degree(), " deg F4'(1,t) =", F4t.degree())
print("Res_t(F3(1,t),F4'(1,t)) =", sp.resultant(F3t, F4t))
print("at direction (0,1): F3 =", F3.subs({x: 0, y: 1}), " F4' =", F4p.subs({x: 0, y: 1}))
c3 = sp.Poly(8*F3.subs({x: 1, y: t}), t)
print("8*F3(1,t) =", c3.as_expr(), " irreducible over Q:", c3.is_irreducible, " disc =", sp.discriminant(c3), sp.factorint(sp.discriminant(c3)))
print("real roots:", sp.Poly(c3).count_roots())
# 4. lines through p tangent to B: discriminant of g_l(m) as quartic in m for direction (1,t)
mm = sp.symbols('mm')
gl = sp.Poly(sp.expand(f.subs({x: mm, y: mm*t, z: 1})), mm)
dl = sp.factor(sp.discriminant(gl, mm))
print("disc of f on line (m, m t, 1):", sp.factor_list(sp.discriminant(gl, mm)))
