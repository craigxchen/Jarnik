# Killer: calibration surface S_0 : w^2 = z^4 + z f3(x,y) + x^4, P0=(0:0:1:1), cuspidal cubic C_0.
import sympy as sp, math
from fractions import Fraction as Fr
x, y, z, e, t = sp.symbols('x y z e t')
f3 = 3*x**3 + 2*x**2*y + x*y**2 + 2*y**3
f = z**4 + z*f3 + x**4
G = sp.groebner([sp.diff(f, v) for v in (x, y, z)], x, y, z, order='grevlex')
lms = [sp.Poly(g, x, y, z).monoms(order='grevlex')[0] for g in G.exprs]
print("pure powers among LMs:", sorted(m for m in lms if sum(1 for u in m if u) == 1))
g = sp.expand(f.subs(z, 1).subs({x: e*x, y: e*y}))
w = sp.expand(sp.series(sp.sqrt(g), e, 0, 5).removeO())
print("local w at P0:", w)
T = z**2  # 2-jet
rem = sp.expand(f - T**2); F3 = rem.coeff(z, 1); F4 = rem.coeff(z, 0)
print("F3 =", F3, " F4' =", F4, " Res =", sp.resultant(sp.Poly(F3.subs({x: 1, y: t}), t), sp.Poly(F4.subs({x: 1, y: t}), t)),
      " F3(0,1) =", F3.subs({x: 0, y: 1}), "F4'(0,1) =", F4.subs({x: 0, y: 1}))
print("disc F3(1,t):", sp.discriminant(sp.Poly(F3.subs({x: 1, y: t}), t)), " irreducible:", sp.Poly(F3.subs({x: 1, y: t}), t).is_irreducible)
# C_0 lies on S_0
print("on C_0: (z^2+x^2)^2 - f restricted to f3=2x^2z:", sp.expand(((z**2 + x**2)**2 - f).subs(z, f3/(2*x**2))))
for bb in (11, 101, 1001, 10001, 100001):
    N = (3 + 2*bb + bb**2 + 2*bb**3)//2
    assert 2*N == 3 + 2*bb + bb**2 + 2*bb**3
    sv = N*N + 1
    assert sv*sv == N**4 + N*(3 + 2*bb + bb*bb + 2*bb**3) + 1
    X = Fr(bb, N); u = sv - N*N
    print(f"b={bb}: c={N}, u={u}, X^1.5*H = {float(X)**1.5*N:.8f}, gcd(a,b)={math.gcd(1,bb)}")
