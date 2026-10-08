"""Exact algebra and height checks for five_row_singular_point_height.md.

Finite fixtures check the algebraic certificates, not a lattice-arc family.
The irreducibility and uniform height argument is proved in the note.
"""

from fractions import Fraction
from itertools import product
from math import gcd

from check_five_row_del_pezzo_arithmetic import graphs
from check_smaller_cut_relation_rank import add, derivative, mul, scale, variable


ONE = {(0, 0): 1}
A, B = variable(2, 0), variable(2, 1)


def evaluate(poly, a, b):
    return sum(value * a**i * b**j for (i, j), value in poly.items())


def l1(poly):
    return sum(abs(value) for value in poly.values())


def coefficient(poly, axis, exponent):
    return {tuple(0 if j == axis else power for j, power in enumerate(mon)): value
            for mon, value in poly.items() if mon[axis] == exponent}


def discriminant(poly, axis):
    quadratic, linear, constant = (coefficient(poly, axis, power)
                                   for power in (2, 1, 0))
    result = add(mul(linear, linear), scale(mul(quadratic, constant), -4))
    identity = add(mul(derivative(poly, axis), derivative(poly, axis)),
                   scale(mul(quadratic, poly), -4))
    assert result == identity
    assert l1(result) <= l1(poly)**2
    return result


def graph_coordinates():
    rows = [({}, ONE), (ONE, {}), (ONE, ONE), (ONE, A), (ONE, B)]
    coordinates = []
    for graph in graphs():
        poly = ONE
        for i, j in graph:
            xi, yi = rows[i]
            xj, yj = rows[j]
            bracket = add(mul(yi, xj), scale(mul(xi, yj), -1))
            poly = mul(poly, bracket)
        coordinates.append(poly)
    return coordinates


def primitive_integer_polynomial(poly):
    denominator = 1
    for value in poly.values():
        denominator = denominator * value.denominator // gcd(denominator, value.denominator)
    integer = {mon: int(value * denominator) for mon, value in poly.items()}
    content = 0
    for value in integer.values():
        content = gcd(content, abs(value))
    return {mon: value // content for mon, value in integer.items() if value}


def main():
    coordinates = graph_coordinates()
    assert len(coordinates) == 22
    for poly in coordinates:
        assert poly and all(i <= 2 and j <= 2 for i, j in poly)
        assert l1(poly) <= 32
        assert not evaluate(poly, 0, 0) and not evaluate(poly, 1, 1)
        assert not poly.get((2, 2), 0)

    # Every coefficient in {-1,0,1} in a five-term moving quadratic.
    monomials = ((2, 0), (1, 2), (1, 1), (0, 2), (0, 0))
    checks = 0
    for values in product((-1, 0, 1), repeat=len(monomials)):
        poly = {mon: value for mon, value in zip(monomials, values) if value}
        if not poly:
            continue
        discriminant(poly, 0)
        discriminant(poly, 1)
        checks += 1

    # A rational nodal fixture of bidegree (2,2), translated to rational
    # singular coordinates. It checks the elementary chart lemma and
    # does not purport to be an anticanonical section in the fixed chart.
    fixtures = 0
    for a0, b0 in ((Fraction(2), Fraction(3)),
                   (Fraction(17, 5), Fraction(-23, 7)),
                   (Fraction(101, 13), Fraction(211, 19))):
        x = add(A, scale(ONE, -a0))
        y = add(B, scale(ONE, -b0))
        poly = primitive_integer_polynomial(add(mul(x, x),
                                                scale(mul(mul(y, y), add(ONE, x)), -1)))
        assert evaluate(poly, a0, b0) == 0
        assert evaluate(derivative(poly, 0), a0, b0) == 0
        assert evaluate(derivative(poly, 1), a0, b0) == 0
        for axis, point in ((0, b0), (1, a0)):
            disc = discriminant(poly, axis)
            assert disc and evaluate(disc, a0, b0) == 0
            assert max(abs(point.numerator), point.denominator) <= l1(disc)
            assert l1(disc) <= l1(poly)**2
        fixtures += 1

    # Exact denominator clearing and the claimed graph-height upper bound.
    chart_fixtures = 0
    rational_values = (Fraction(-7, 3), Fraction(2), Fraction(3),
                       Fraction(17, 5), Fraction(29, 11))
    for a, b in product(rational_values, repeat=2):
        if a == b:
            continue
        denominator = a.denominator**2 * b.denominator**2
        values = [evaluate(poly, a, b) * denominator for poly in coordinates]
        assert all(value.denominator == 1 and value for value in values)
        bound = (32 * max(abs(a.numerator), a.denominator)**2
                 * max(abs(b.numerator), b.denominator)**2)
        assert max(abs(value) for value in values) <= bound
        chart_fixtures += 1

    # A smooth rational ramification point on the known anticanonical
    # elliptic section G_12; it is not a singular-point fixture.
    elliptic = {(2, 0): -12, (1, 0): 12, (1, 1): 13,
                (1, 2): -1, (0, 1): -12}
    branch = discriminant(elliptic, 0)
    expected = mul(mul(add(B, scale(ONE, -1)), add(B, scale(ONE, -4))),
                   add(add(mul(B, B), scale(B, -21)), scale(ONE, 36)))
    assert branch == expected
    assert evaluate(elliptic, 2, 4) == evaluate(derivative(elliptic, 0), 2, 4) == 0
    assert evaluate(derivative(elliptic, 1), 2, 4) == -2
    leading = coefficient(elliptic, 0, 2)
    linear = coefficient(elliptic, 0, 1)
    recovered = -Fraction(evaluate(linear, 2, 4), 2 * evaluate(leading, 2, 4))
    assert recovered == 2
    assert max(abs(recovered.numerator), recovered.denominator) <= 2*l1(elliptic)*4**2

    print(f"PASS: 22 graph charts; bidegrees, section constraints, l1 max {max(map(l1, coordinates))}.")
    print(f"PASS: {checks} exact discriminant identities in both variables; coefficient-square bounds.")
    print(f"PASS: {fixtures} rational singular fixtures and {chart_fixtures} exact chart-height checks.")
    print("PASS: smooth rational ramification fixture, discriminant and coordinate recovery.")
    print("No endpoint family or exclusion of unramified smooth points is asserted.")


if __name__ == "__main__":
    main()
