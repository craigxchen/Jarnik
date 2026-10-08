"""Exact five-row normalized graph products and an interior nodal section."""

from fractions import Fraction
from math import gcd

from check_five_row_del_pezzo_arithmetic import graphs
from check_five_row_boundary_slope_relations import EXPECTED


ONE = {(0, 0): 1}
A = {(1, 0): 1}
B = {(0, 1): 1}


def add(left, right):
    out = dict(left)
    for monomial, coefficient in right.items():
        out[monomial] = out.get(monomial, 0) + coefficient
    return {m: c for m, c in out.items() if c}


def scale(poly, coefficient):
    return {m: coefficient * c for m, c in poly.items() if coefficient * c}


def mul(left, right):
    out = {}
    for (i, j), left_coefficient in left.items():
        for (k, ell), right_coefficient in right.items():
            monomial = (i + k, j + ell)
            out[monomial] = out.get(monomial, 0) + left_coefficient * right_coefficient
    return {m: c for m, c in out.items() if c}


def derivative(poly, variable):
    out = {}
    for (i, j), coefficient in poly.items():
        exponent = (i, j)[variable]
        if exponent:
            monomial = (i - 1, j) if variable == 0 else (i, j - 1)
            out[monomial] = coefficient * exponent
    return out


def evaluate(poly, a, b):
    return sum(coefficient * a**i * b**j
               for (i, j), coefficient in poly.items())


def degree(poly):
    return (max(i for i, _ in poly), max(j for _, j in poly))


def primitive_content(poly):
    content = 0
    for coefficient in poly.values():
        content = gcd(content, abs(coefficient))
    return content


def discriminant_in_a(poly):
    coefficients = [
        {(0, j): value for (i, j), value in poly.items() if i == power}
        for power in (2, 1, 0)
    ]
    quadratic, linear, constant = coefficients
    disc = add(mul(linear, linear), scale(mul(quadratic, constant), -4))
    return {j: value for (_, j), value in disc.items()}


def discriminant_is_not_square(poly):
    # F=(24+12b)a^2+(-3-58b-7b^2)a+(8b+24b^2).
    # Its discriminant is 49b^4-340b^3+718b^2-420b+9.
    discriminant = {0: 9, 1: -420, 2: 718, 3: -340, 4: 49}
    assert discriminant_in_a(poly) == discriminant
    quadratic_factor = {0: 1, 1: -46, 2: 49}
    squared_factor = {0: 9, 1: -6, 2: 1}
    factored = {}
    for degree, coefficient in squared_factor.items():
        for other_degree, other_coefficient in quadratic_factor.items():
            factored[degree + other_degree] = (
                factored.get(degree + other_degree, 0)
                + coefficient * other_coefficient
            )
    assert factored == discriminant
    assert (-46) ** 2 - 4 * 49 * 1 == 1920
    assert discriminant[4] == 49 and discriminant[0] == 9
    # The factored quadratic has two simple roots over Qbar, so the quartic
    # has simple roots away from b=3 and cannot be a square in Qbar(b).


def main():
    # Rows are (0,1),(1,0),(1,1),(1,a),(1,b), with
    # Delta_ij=Y_i X_j-X_i Y_j.
    brackets = {
        (0, 1): ONE, (0, 2): ONE, (0, 3): ONE, (0, 4): ONE,
        (1, 2): scale(ONE, -1), (1, 3): scale(A, -1),
        (1, 4): scale(B, -1), (2, 3): add(ONE, scale(A, -1)),
        (2, 4): add(ONE, scale(B, -1)),
    }
    brackets[(3, 4)] = add(A, scale(B, -1))
    products = []
    for graph in graphs():
        product = ONE
        for i, j in graph:
            key = tuple(sorted((i, j)))
            factor = brackets[key] if i < j else scale(brackets[key], -1)
            product = mul(product, factor)
        products.append(product)

    assert len(products) == 22
    assert all(degree(product)[0] <= 2 and degree(product)[1] <= 2
               for product in products)
    assert all(product.get((2, 2), 0) == 0 for product in products)
    assert all(evaluate(product, 0, 0) == 0 for product in products)
    assert all(evaluate(product, 1, 1) == 0 for product in products)
    l1 = [sum(abs(c) for c in product.values()) for product in products]
    coefficient_max = [max(abs(c) for c in product.values()) for product in products]
    assert max(l1) == 6 and max(coefficient_max) == 2
    assert max(l1) <= 32 and max(coefficient_max) <= 32

    # Primitive section with an interior singular point (a,b)=(2,3).
    section_coefficients = (-24, 7, -12, 0, 0, 4)
    section = {}
    for coefficient, product in zip(section_coefficients, products[:6]):
        section = add(section, scale(product, coefficient))
    assert section == {
        (2, 0): -24, (1, 1): 58, (0, 2): -24,
        (1, 0): 3, (1, 2): 7, (0, 1): -8, (2, 1): -12,
    }
    assert primitive_content(section) == 1
    assert evaluate(section, 2, 3) == 0
    assert evaluate(derivative(section, 0), 2, 3) == 0
    assert evaluate(derivative(section, 1), 2, 3) == 0
    assert discriminant_is_not_square(section) is None

    # Gauss primitivity over Qbar[b] is separate from integer coefficient
    # content.  The leading and linear-in-a coefficients are
    # A=12(b+2), B=-3-58b-7b^2, and B(-2)=85, so gcd(A,B)=1.
    assert -3 - 58 * (-2) - 7 * (-2) ** 2 == 85

    # The Hessian is nonsingular, so this singularity is an ordinary node.
    hessian = (
        (evaluate(derivative(derivative(section, 0), 0), 2, 3),
         evaluate(derivative(derivative(section, 0), 1), 2, 3)),
        (evaluate(derivative(derivative(section, 1), 0), 2, 3),
         evaluate(derivative(derivative(section, 1), 1), 2, 3)),
    )
    assert hessian == ((-120, 52), (52, -20))
    assert hessian[0][0] * hessian[1][1] - hessian[0][1] * hessian[1][0] == -304

    # A second primitive singular section has no nodal boundary contact.
    # This is the vector (-1,-1,0) in the three-dimensional singular kernel.
    extra_coefficients = (116, -33, 72, -4, -4, 0)
    extra = {}
    for coefficient, product in zip(extra_coefficients, products[:6]):
        extra = add(extra, scale(product, coefficient))
    assert extra == {
        (2, 0): 120, (1, 1): -310, (0, 2): 120,
        (1, 0): -33, (1, 2): -37, (0, 1): 72, (2, 1): 68,
    }
    assert primitive_content(extra) == 1
    assert evaluate(extra, 2, 3) == 0
    assert evaluate(derivative(extra, 0), 2, 3) == 0
    assert evaluate(derivative(extra, 1), 2, 3) == 0
    extra_hessian = (
        (evaluate(derivative(derivative(extra, 0), 0), 2, 3),
         evaluate(derivative(derivative(extra, 0), 1), 2, 3)),
        (evaluate(derivative(derivative(extra, 1), 0), 2, 3),
         evaluate(derivative(derivative(extra, 1), 1), 2, 3)),
    )
    assert extra_hessian == ((648, -260), (-260, 92))
    assert extra_hessian[0][0] * extra_hessian[1][1] - extra_hessian[0][1] * extra_hessian[1][0] == -7984
    for alpha, beta in EXPECTED.values():
        av = sum(coefficient * value for coefficient, value in zip(alpha, extra_coefficients))
        bv = sum(coefficient * value for coefficient, value in zip(beta, extra_coefficients))
        assert av and bv and av + bv

    # Its a-discriminant factors as (b-3)^2 times a quadratic with distinct
    # roots over Qbar, so it is nonsquare in Qbar(b).  Primitivity follows
    # from A=68b+120 and B(-30/17)=115263/289 != 0.
    extra_discriminant = {0: 1089, 1: -14100, 2: 21358, 3: -9700, 4: 1369}
    assert discriminant_in_a(extra) == extra_discriminant
    extra_factor = {0: 121, 1: -1486, 2: 1369}
    assert extra_discriminant == {
        degree: sum((9 if degree - shift == 0 else -6 if degree - shift == 1 else 1)
                    * value for shift, value in extra_factor.items()
                    if 0 <= degree - shift <= 2)
        for degree in range(5)
    }
    assert (-1486) ** 2 - 4 * 1369 * 121 == 1545600
    root = Fraction(-30, 17)
    assert -33 - 310 * root - 37 * root**2 == Fraction(115263, 289)

    # The three prescribed blown-up points have multiplicity exactly one;
    # thus homogenization and blow-up introduce no exceptional component.
    for point in ((0, 0), (1, 1)):
        assert evaluate(extra, *point) == 0
        assert any(evaluate(derivative(extra, axis), *point) for axis in (0, 1))
    at_infinity = {(2-i, 2-j): value for (i, j), value in extra.items()}
    assert evaluate(at_infinity, 0, 0) == 0
    assert any(evaluate(derivative(at_infinity, axis), 0, 0) for axis in (0, 1))

    # (2,3) avoids all boundary values 0,1,infinity and a=b.
    assert 2 not in (0, 1) and 3 not in (0, 1) and 2 != 3
    print("PASS: all 22 normalized graph products have bidegree <=(2,2),")
    print("      vanish at (0,0),(1,1), and have zero a^2 b^2 coefficient.")
    print(f"      max l1={max(l1)}, max coefficient={max(coefficient_max)} (bounds 32,32).")
    print("PASS: primitive singular F fixtures include a nodal one with all ten boundary restrictions nonnodal.")
    print("No uniform arc or family theorem is asserted.")


if __name__ == "__main__":
    main()
