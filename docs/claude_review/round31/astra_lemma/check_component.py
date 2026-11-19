"""Exact identities and the 2-adic obstruction for the round-24 family M."""

import sympy as sp


m, y = sp.symbols('m y')
c = m**3 - 11*m**2 + 39*m - 13
D = m**4 + 70*m**2 - 48*m + 9
A = -(m + 3)*c
B = 2*D
u = A / (16*(m - 3)**2*(m - 1)**2)
v = -(m + 1)*c / (4*(m - 3)**3*(m - 1))
G = (y**2 - y + u)**2 + v*y

# Dividing the cleared quartic discriminant by its leading coefficient to
# power six recovers the discriminant of the monic quartic G.
num, den = sp.together(G).as_numer_denom()
disc = sp.discriminant(num, y) / den**6
factor = 2*(m + 1)*c / ((m - 3)**5*(m - 1)**3)
assert sp.cancel(disc - u*B*factor**2) == 0
assert sp.cancel(v/u - 4*(m**2 - 1)/(m**2 - 9)) == 0
assert sp.cancel(u - u.subs(m, -m) + 256*m**3/((m**2 - 9)**2*(m**2 - 1)**2)) == 0
assert sp.limit(u, m, sp.oo) == -sp.Rational(1, 16)

a = sp.lambdify(m, A)
b = sp.lambdify(m, B)
squares = {j*j % 256 for j in range(256)}
for r in (1, 3, 5, 7):
    values = {(a(j) % 256, b(j) % 256) for j in range(r, 256, 8)}
    assert len(values) == 1
    first, second = values.pop()
    assert first not in squares or second not in squares
    print(f'm = {r} mod 8: (A,B) = ({first},{second}) mod 256')
print('Quartic discriminant, parameter recovery and all odd 2-adic residues verified.')

# Split-form control: both square classes change sign. The full splitting
# condition at 3 then excludes the parameter chart at infinity.
squares_81 = {j*j % 81 for j in range(81)}
assert all(-a(j) % 81 not in squares_81 or -b(j) % 81 not in squares_81
           for j in range(81))
g_at_infinity = sp.Poly(sp.limit(G, m, sp.oo), y)
coefficients = [int(c.p) * pow(int(c.q), -1, 3) % 3
                for c in g_at_infinity.all_coeffs()]
assert coefficients == [1, 1, 2, 1, 1]
g_reduced = sp.Poly.from_list(coefficients, y, modulus=3)
assert [j for j in range(3) if g_reduced.eval(j) % 3 == 0] == [1]
assert g_reduced != sp.Poly((y - 1)**4, y, modulus=3)
print('Split-form control at 3 verified.')
