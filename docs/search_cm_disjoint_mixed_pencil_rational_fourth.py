"""Optional SymPy search: three integral common roots and a rational fourth."""

from itertools import combinations

import sympy as s


X = s.symbols("X")
K = s.QQ.frac_field(X)


def branch_functions(triple):
    a, b, c = map(K, triple)
    x = K.gens[0]
    total = a + b + c
    pairs = a*b + a*c + b*c
    product = a*b*c
    p0 = product*x
    p1 = -product + x*(product - pairs)
    q1 = -total - x
    q0 = pairs - product + x*(total + product - pairs)
    lambda_b = 4*p0/p1**2
    u, v = p1 - 4*q1, p0 - 4*q0
    lambda_c = (8*u + 4*v - 64)/(u**2 + 16*v)

    def discriminant_coefficients(lam):
        d1 = (-2*lam*p1*(1 + lam*q1) + 4*(-1 + lam*q0)
              + 4*lam*lam*p0)
        d2 = (1 + lam*q1)**2 - 4*lam*(-1 + lam*q0)
        return d1, d2

    d1, d2 = discriminant_coefficients(lambda_b)
    branch_s = -d1/d2
    d1, d2 = discriminant_coefficients(lambda_c)
    branch_t = -d1/d2 - 4
    return branch_s, branch_t


def candidates(bound=15):
    values = [n for n in range(-bound, bound + 1) if n not in (0, 1, 2)]
    for triple in combinations(values, 3):
        try:
            branch_s, branch_t = branch_functions(triple)
        except ZeroDivisionError:
            continue
        harmonic = (branch_s*branch_t - 2*branch_s - 2*branch_t,
                    branch_s*branch_t + 4*branch_s - 8*branch_t,
                    branch_s*branch_t - 8*branch_s + 4*branch_t)
        for ordering, equation in enumerate(harmonic):
            numerator = s.Poly(equation.numer.as_expr(), X)
            for fourth in s.polys.polytools.ground_roots(numerator):
                if fourth.is_Rational and fourth not in (*triple, 0, 1, 2):
                    yield triple, fourth, ordering


if __name__ == "__main__":
    found = list(candidates())
    assert not found, found
    print("No harmonic rational fourth root occurs for three integral roots in the box.")
