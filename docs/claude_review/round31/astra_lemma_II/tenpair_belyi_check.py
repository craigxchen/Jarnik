"""Exact residual critical-value calculation for the plane-pencil proof."""
import sympy as sp

x, b = sp.symbols("x b")
q = 2*x**2 + x - 3*b
num, den = x**3, (x - 1)**3*(x - b)**2
critical = sp.diff(num, x)*den - num*sp.diff(den, x)
assert sp.expand(critical + x**2*(x-1)**2*(x-b)*q) == 0
a0, a1 = [sp.rem(num, q, x).expand().coeff(x, j) for j in (0, 1)]
c0, c1 = [sp.rem(den, q, x).expand().coeff(x, j) for j in (0, 1)]
det = sp.factor(a0*c1-a1*c0)
assert sp.factor(16*det - b**2*(9*b+16)*(24*b+1)) == 0
assert sp.discriminant(q, x) == 24*b + 1
print("constant-linear determinant =", det)
print("discriminant at -16/9 =", (24*b+1).subs(b, -sp.Rational(16,9)))
