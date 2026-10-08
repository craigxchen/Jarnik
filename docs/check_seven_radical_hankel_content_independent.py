#!/usr/bin/env python3
"""Independent exact checks for the seven-row radical Hankel identities."""

from __future__ import annotations

from fractions import Fraction as F
from functools import reduce
from itertools import combinations, product
from math import gcd, isqrt, lcm, prod
from random import Random


LABELS = range(7)
PAIRS = list(combinations(LABELS, 2))


def determinant(matrix: list[list[F]]) -> F:
    work = [row[:] for row in matrix]
    answer = F(1)
    for column in range(len(work)):
        pivot = next(row for row in range(column, len(work)) if work[row][column])
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            answer = -answer
        value = work[column][column]
        answer *= value
        work[column] = [entry / value for entry in work[column]]
        for row in range(column + 1, len(work)):
            multiple = work[row][column]
            work[row] = [x - multiple * y for x, y in zip(work[row], work[column])]
    return answer


def sqrt_fraction(value: F) -> F:
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    assert numerator * numerator == value.numerator
    assert denominator * denominator == value.denominator
    return F(numerator, denominator)


def rational_gcd(values: list[F]) -> F:
    return F(
        reduce(gcd, (abs(value.numerator) for value in values)),
        reduce(lcm, (value.denominator for value in values), 1),
    )


def affine_matrix(ts: list[F], ys: list[F]) -> list[list[F]]:
    pprime = [prod(ts[i] - ts[j] for j in LABELS if i != j) for i in LABELS]
    moments = [sum(ys[i] * ts[i] ** degree / pprime[i] for i in LABELS)
               for degree in range(5)]
    return [[moments[a + b] for b in range(3)] for a in range(3)]


def check_affine_branches() -> None:
    tuples = [
        [F(j) for j in range(1, 8)],
        [F(j, 2) for j in range(1, 8)],
        [F(j * j + 1, j + 1) for j in range(1, 8)],
    ]
    for initial_us in tuples:
        ts = [(u - 1 / u) / 2 for u in initial_us]
        assert len(set(ts)) == 7
        vandermonde = prod(ts[j] - ts[i] for i, j in PAIRS)
        branch_product = F(1)
        for signs in product((-1, 1), repeat=7):
            ys0 = [(u + 1 / u) / 2 for u in initial_us]
            ys = [signs[i] * ys0[i] for i in LABELS]
            us = [ys[i] + ts[i] for i in LABELS]
            actual = determinant(affine_matrix(ts, ys))
            expected = -F(2**12) * prod(u**3 for u in us) / prod(
                1 + us[i] * us[j] for i, j in PAIRS
            )
            assert actual == expected and actual
            pair_denominator = prod(
                1 + ts[i] * ts[j] + ys[i] * ys[j] for i, j in PAIRS
            )
            assert actual * actual == F(8, pair_denominator)
            branch_product *= actual
        assert branch_product == F(2**192, vandermonde**64)


def configuration(form: tuple[F, F, F], vectors: list[tuple[int, int]],
                  check_radicals: bool = False) -> dict[str, object]:
    a, b, c = form
    radius = sqrt_fraction(a * c - b * b)
    qs = [a*x*x + 2*b*x*y + c*y*y for x, y in vectors]
    brackets = {(i, j): F(vectors[i][0] * vectors[j][1]
                                - vectors[i][1] * vectors[j][0])
                for i, j in PAIRS}

    def bracket(i: int, j: int) -> F:
        return brackets[i, j] if i < j else -brackets[j, i]

    stars = [prod(bracket(i, j) for j in LABELS if j != i) for i in LABELS]
    D = abs(prod(brackets.values()))
    rho_squares = [qs[i]**5 / stars[i]**2 for i in LABELS]
    tau = rational_gcd(rho_squares)
    cs = [value / tau for value in rho_squares]
    assert all(value.denominator == 1 and value > 0 for value in cs)

    triple_values = []
    term_squares = []
    for subset in combinations(LABELS, 3):
        complement = [i for i in LABELS if i not in subset]
        vi = prod(brackets[i, j] for i, j in combinations(subset, 2))
        vc = prod(brackets[i, j] for i, j in combinations(complement, 2))
        triple_values.append(prod(qs[i] for i in subset) * vi**2 * vc**2)
        term_squares.append(
            prod(cs[i] for i in subset) * vi**4 / prod(qs[i]**4 for i in subset)
        )
    G = rational_gcd(triple_values)
    g3 = rational_gcd(term_squares)
    assert g3 == G / (tau**3 * D**2)
    base = 8 * radius**3 * D / G
    assert base.denominator == 1 and base > 0

    result: dict[str, object] = {"base": base, "cs": cs, "tau": tau}
    if not check_radicals:
        return result

    roots = [sqrt_fraction(value) for value in cs]
    signs = [1 if star > 0 else -1 for star in stars]
    features = [[F(x*x), F(x*y), F(y*y)] for x, y in vectors]
    matrix = [[sum(signs[i] * roots[i] * features[i][j] * features[i][k] / qs[i]**2
                          for i in LABELS)
               for k in range(3)] for j in range(3)]
    z3 = determinant(matrix) / sqrt_fraction(g3)
    assert z3.denominator == 1

    unit_product = F(1)
    for i, j in PAIRS:
        x, y = vectors[i]
        z, t = vectors[j]
        bilinear = a*x*z + b*(x*t + y*z) + c*y*t
        root_product = sqrt_fraction(qs[i] * qs[j])
        pair = (bilinear + root_product) / (radius * abs(brackets[i, j]))
        conjugate = (bilinear - root_product) / (radius * abs(brackets[i, j]))
        assert pair > 0 and pair * conjugate == -1
        unit_product *= pair
    assert z3*z3 == base / unit_product

    formal_product = F(1)
    for flips in product((-1, 1), repeat=7):
        signed_matrix = [[
            sum(flips[i] * signs[i] * roots[i]
                * features[i][j] * features[i][k] / qs[i]**2 for i in LABELS)
            for k in range(3)] for j in range(3)]
        signed_z3 = determinant(signed_matrix) / sqrt_fraction(g3)
        assert signed_z3.denominator == 1
        formal_product *= signed_z3
    assert formal_product == base**64
    result.update({"z3": z3, "units": unit_product})
    return result


def check_content_covariance() -> None:
    rng = Random(20260913)
    forms = [(F(1), F(0), F(1)), (F(5), F(2), F(1)),
             (F(2), F(0), F(2))]
    checked = 0
    for form in forms:
        for _ in range(40):
            while True:
                vectors = [(rng.randrange(1, 10), rng.randrange(-12, 13))
                           for _ in LABELS]
                if all(vectors[i][0] * vectors[j][1] !=
                       vectors[i][1] * vectors[j][0] for i, j in PAIRS):
                    break
            configuration(form, vectors)
            checked += 1

    # Rational conic parametrization makes every required radical rational.
    vectors = [(2*u, u*u - 1) for u in range(1, 8)]
    original = configuration((F(1), F(0), F(1)), vectors, True)
    diagonal_vectors = [(2*x, 3*y) for x, y in vectors]
    diagonal = configuration((F(1, 4), F(0), F(1, 9)), diagonal_vectors, True)
    shear_vectors = [(x + 2*y, y) for x, y in vectors]
    shear = configuration((F(1), F(-2), F(5)), shear_vectors, True)
    scaled = configuration((F(2, 3), F(0), F(2, 3)), vectors, True)
    assert original["base"] == diagonal["base"] == shear["base"] == scaled["base"]
    assert original["cs"] == diagonal["cs"] == shear["cs"] == scaled["cs"]
    assert original["z3"] == diagonal["z3"] == shear["z3"] == scaled["z3"]
    assert original["units"] == diagonal["units"] == shear["units"] == scaled["units"]
    print(f"checked integer base on {checked} random marked configurations")


def poly_multiply(left: list[F], right: list[F]) -> list[F]:
    answer = [F(0)] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            answer[i + j] += x * y
    return answer


def check_individual_forms() -> None:
    vectors = [(1, t) for t in (-7, -4, -2, 1, 3, 6, 10)]
    form = (F(5), F(2), F(1))
    result = configuration(form, vectors)
    cs = result["cs"]
    a, b, c = form
    qs = [a*x*x + 2*b*x*y + c*y*y for x, y in vectors]
    for r in range(4):
        terms = []
        for i, (x, y) in enumerate(vectors):
            for j in range(2*r + 1):
                terms.append(cs[i] * F(x)**(2*(2*r-j)) * F(y)**(2*j)
                             / qs[i]**(2*r))
        predicted = rational_gcd([cs[i] / qs[i]**(2*r) for i in LABELS])
        assert rational_gcd(terms) == predicted

    ts = [F(y, x) for x, y in vectors]
    source = [F(i*i + 3*i + 1, i + 2) for i in LABELS]
    forms = [sum(source[i] * ts[i]**j / qs[i]**3 for i in LABELS)
             for j in range(7)]
    for i in LABELS:
        numerator = [F(1)]
        denominator = F(1)
        for j in LABELS:
            if i != j:
                numerator = poly_multiply(numerator, [-ts[j], F(1)])
                denominator *= ts[i] - ts[j]
        coefficients = [value / denominator for value in numerator]
        assert coefficients[6] == 1 / denominator
        assert qs[i]**3 * sum(coefficients[j] * forms[j] for j in range(7)) == source[i]


def primitive_gaussian(real: F, imag: F) -> tuple[int, int]:
    denominator = lcm(real.denominator, imag.denominator)
    x = int(real * denominator)
    y = int(imag * denominator)
    content = gcd(abs(x), abs(y))
    return x // content, y // content


def check_squareclasses() -> None:
    form = (F(5), F(2), F(1))  # determinant one
    vectors = [(1, -4), (2, -3), (3, -1), (4, 1), (3, 2), (2, 3), (1, 5)]
    result = configuration(form, vectors)
    cs = result["cs"]
    tau = result["tau"]
    a, b, c = form
    r = sqrt_fraction(a*c - b*b)
    qs = [a*x*x + 2*b*x*y + c*y*y for x, y in vectors]
    brackets = {(i, j): F(vectors[i][0] * vectors[j][1]
                                - vectors[i][1] * vectors[j][0])
                for i, j in PAIRS}
    stars = []
    for i in LABELS:
        stars.append(prod(brackets[i, j] if i < j else -brackets[j, i]
                          for j in LABELS if i != j))

    ells = [(a*x + b*y, r*y) for x, y in vectors]
    for i, j in PAIRS:
        # beta=ell_i conjugate(ell_j), then clear its arbitrary Q-content.
        real = ells[i][0] * ells[j][0] + ells[i][1] * ells[j][1]
        imag = ells[i][1] * ells[j][0] - ells[i][0] * ells[j][1]
        hx, hy = primitive_gaussian(real, imag)
        norm_h = F(hx*hx + hy*hy)
        sqrt_fraction(norm_h / (qs[i] * qs[j]))
        sqrt_fraction(cs[i] * cs[j] / (qs[i] * qs[j]))
        expected_square = qs[i]**2 * qs[j]**2 / (stars[i] * stars[j] * tau)
        assert cs[i] * cs[j] / (qs[i] * qs[j]) == expected_square**2

    # (1+i)/(1-i)=i, so the full unit ratio carries norm class 2.
    numerator, conjugate = (1, 1), (1, -1)
    i_times_conjugate = (-conjugate[1], conjugate[0])
    assert numerator == i_times_conjugate
    assert numerator[0]**2 + numerator[1]**2 == 2


def main() -> None:
    check_affine_branches()
    check_content_covariance()
    check_individual_forms()
    check_squareclasses()
    print("PASS: affine constants and all 128 formal branches.")
    print("PASS: complete contents, integer base, pair factors, and covariance.")
    print("PASS: full half-angle norm and primitive Gale squareclasses agree pairwise.")


if __name__ == "__main__":
    main()
