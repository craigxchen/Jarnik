"""(1) Universal weighted family: s^2 = c^4 + 2N c m^3 + m^4 has the point (m,c,s) = (1,N,N^2+1); weighted height
A' = (2N)^{1/3}, X = A'/N ~ c^{-2/3}: Liouville exponent 3/2 attained, uniformly in N.
(2) Its realisation on a special dP2: S0: w^2 = z^4 + z*f3(x,y) + x^4, P0 = (0:0:1:1). The points
Q_b = (1 : b : N : N^2+1), N = f3(1,b)/2, lie on the cuspidal cubic C0 = {f3(x,y) = 2x^2 z, s = z^2 + x^2},
with u = s - T(Q) = 1 for all b. Check smoothness of S0, P0 off lines (Res(F3,F4') != 0), and X^{3/2} H bounded."""
import sympy as sp
from sympy import symbols, expand, Poly, groebner, diff, resultant
from geom import analyse
x, y, z, t = symbols('x y z t')
f3 = 3*x**3 + 2*x**2*y + x*y**2 + 2*y**3      # f3(1,b) = 3 + 2b + b^2 + 2b^3: even iff b odd
f0 = expand(z**4 + z*f3 + x**4)
G = groebner([diff(f0, v) for v in (x, y, z)], x, y, z, order='grevlex')
lms = [Poly(g, x, y, z).monoms(order='grevlex')[0] for g in G.exprs]
pure = {i for mon in lms for i in range(3) if mon[i] > 0 and sum(1 for e in mon if e > 0) == 1}
print('S0 smooth over Q:', pure == {0, 1, 2})
out, w, T2hat, F4 = analyse(f0)
print('S0: T =', T2hat, '| F3 =', out['F3'], '| F4prime =', out['F4prime'], '| Res(F3,F4p) =', out['Res_F3_F4p'],
      '| F3 factor:', out['F3_factor'])
for b in [3, 11, 101, 1001, 10001, 100001]:
    N = f3.subs({x: 1, y: b}) / 2
    assert N == int(N)
    N = int(N)
    a_, b_, c_ = 1, b, N
    s_ = N**2 + 1
    assert f0.subs({x: a_, y: b_, z: c_}) == s_**2
    u = s_ - T2hat.subs({x: a_, y: b_, z: c_})
    X = max(abs(a_), abs(b_)) / abs(c_)
    print('b=%d: H=c=%d, u=%s, X^{3/2} H = %.6f' % (b, c_, u, X**1.5 * c_))
