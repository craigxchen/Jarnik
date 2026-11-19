"""Rational certificates for the four isotropic star blocks, Part V."""
from itertools import combinations

import sympy as sp


v, t, k = sp.symbols('v t k', nonzero=True)
r = t**2
f = (1-r)*(r**2-3*r+1)
x = (k+f/k)/2
y = (f/k-k)/(2*t)
roots = [v*(1+x)/2, v*(1-x)/2, v*(-r+y)/2, v*(-r-y)/2]
point = {v: sp.Rational(3, 2), t: 2, k: 2}
aa = [a.subs(point) for a in roots]
us = [1] + [a*a for a in aa]
assert len({us[i]*us[j] for i, j in combinations(range(5), 2)}) == 10
s1 = sum(aa)
s2 = sum(a*b for a, b in combinations(aa, 2))
s3 = sum(a*b*c for a, b, c in combinations(aa, 3))
s4 = sp.prod(aa)
assert s1*s2 == s3 and s1*s4 != 0
direction = [a**(-2)-1 for a in roots]
jac = sp.Matrix([[d.subs(point)] + [sp.diff(d, z).subs(point) for z in (v, t, k)]
                 for d in direction])
assert jac.det() != 0
K = sp.symbols('K', real=True)
C = sp.expand(sp.prod(K-sp.I*a for a in aa))
A = K-sp.I*s1
assert sp.expand(sp.im(sp.conjugate(C)*A)) == -s1*s4
W = [u*C+(u-1)*(u+K**2)*A for u in us]
for i, j in combinations(range(5), 2):
    J = K**2+us[i]*us[j]
    assert sp.expand(sp.im(sp.conjugate(W[i])*W[j]) + s1*s4*(us[j]-us[i])*J) == 0
    divis = sp.rem(sp.expand(sp.conjugate(W[i])*W[i]), J, K) == 0
    assert divis == (i == 0)
print('satellites:', aa)
print('projective Jacobian certificate:', jac.det())
