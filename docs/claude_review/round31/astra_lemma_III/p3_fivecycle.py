"""Check the ordinary cyclic order-five symmetry over Q(zeta_5)."""
from sympy import Symbol, cancel, cyclotomic_poly, fraction, invert, rem, resultant

z = Symbol("z")
f = cyclotomic_poly(5, z)

def reduce_rational(expression):
    n, d = fraction(cancel(expression))
    return rem(n*invert(d, f, z), f, z).expand()

G = [0] + [reduce_rational((1+z**((-k)%5))/(1+z**((-2*k)%5)))
           for k in range(1, 5)]
g = G[1]
c = reduce_rational(G[2]/g)
mean_sq = reduce_rational((g**3-G[3])/(3*(G[3]-g*c*c)))
residual = reduce_rational(G[4]*(mean_sq+1)-g*c*(c*c*mean_sq+g*g))
assert mean_sq == 5-8*z**2-8*z**3
assert residual == -24*z**3-32*z**2-24*z
assert resultant(f, mean_sq, z) == 1
assert resultant(f, 3*z**2+4*z+3, z) == 25
print("mean squared", mean_sq, "fourth-moment residual", residual)
