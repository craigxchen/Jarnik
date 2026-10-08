"""Exact supplemental checks for Sections 7--9 of pell_product_endpoint_obstruction.md."""

from fractions import Fraction
from itertools import combinations, product
from math import gcd


def fibonacci(n: int) -> int:
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def convolve(left: list[int], right: list[int]) -> list[int]:
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return out


def weight(coefficients: list[int]) -> int:
    return sum(abs(value) for value in coefficients)


def check_fibonacci_orbit() -> None:
    old_x, old_y = 1, 0
    for r in range(101):
        a, b = fibonacci(2 * r + 1), fibonacci(2 * r)
        assert gcd(a, b) == 1
        assert a * a + b * b == fibonacci(4 * r + 1)
        if r % 3 == 0:
            assert (a, b) == (old_x + old_y, 2 * old_y)
            old_x, old_y = 9 * old_x + 20 * old_y, 4 * old_x + 9 * old_y
        next_a, next_b = fibonacci(2 * r + 3), fibonacci(2 * r + 2)
        assert (next_a, next_b) == (2 * a + b, a + b)

    for n in range(10, 81):
        for s, t in combinations(range(6), 2):
            left = fibonacci(4 * (n + s) + 1)
            right = fibonacci(4 * (n + t) + 1)
            assert fibonacci(4 * (t - s)) % gcd(left, right) == 0


def check_low_height_multiples() -> None:
    minimal_polynomial = [1, -7, 1]
    assert weight(minimal_polynomial) == 9
    assert sum(minimal_polynomial) == -5

    # Translation and overall sign do not affect the coefficient weight.
    # These are all two-unit-monomial relative shifts up to 20.
    smallest_two_term_weight = 100
    for gap in range(1, 21):
        for sign in (-1, 1):
            quotient = [0] * (gap + 1)
            quotient[0], quotient[gap] = 1, sign
            smallest_two_term_weight = min(
                smallest_two_term_weight,
                weight(convolve(minimal_polynomial, quotient)),
            )
    assert smallest_two_term_weight == 14

    # Exhaust all low-weight quotients on a short support; the proof's
    # reverse-triangle estimate covers arbitrary support and coefficients.
    checked = 0
    for coefficients in product(range(-2, 3), repeat=6):
        quotient = list(coefficients)
        quotient_weight = weight(quotient)
        if quotient_weight == 0 or quotient_weight > 3:
            continue
        multiple = convolve(minimal_polynomial, quotient)
        assert weight(multiple) >= 5 * quotient_weight
        if weight(multiple) <= 12:
            assert quotient_weight == 1
            assert sum(multiple) % 2 == 1
        checked += 1
    print(f"PASS: {checked} low-weight quotient fixtures; minimum two-term cost 14")


def check_general_quadratic_units() -> None:
    # (D, Re F, Im F, x_1, y_1, Norm lambda, clearing integer).
    cases = (
        (2, 1, 1, Fraction(3), Fraction(2), 1, 1),
        (2, -1, 1, Fraction(1), Fraction(1), -1, 1),
        (5, 1, 2, Fraction(3, 2), Fraction(1, 2), 1, 2),
        (5, -1, 2, Fraction(3, 2), Fraction(1, 2), 1, 2),
        (5, 1, -2, Fraction(1, 2), Fraction(1, 2), -1, 2),
        (10, 1, 3, Fraction(19), Fraction(6), 1, 1),
        (13, -2, -3, Fraction(649), Fraction(180), 1, 1),
        (13, 2, -3, Fraction(18), Fraction(5), -1, 1),
    )
    fixtures = 0
    for D, a, b, x1, y1, sign, clearing in cases:
        assert a * a + b * b == D
        assert x1 * x1 - D * y1 * y1 == sign
        trace = 2 * x1
        assert trace.denominator == 1
        recurrence_trace = int(trace * trace - 2 * sign)
        xy = [(Fraction(1), Fraction(0))]
        for _ in range(50):
            x, y = xy[-1]
            xy.append((x * x1 + D * y * y1, x * y1 + y * x1))

        norms = []
        for r, (x, y) in enumerate(xy[:25]):
            X, Y = clearing * x, clearing * y
            assert X.denominator == Y.denominator == 1
            X, Y = int(X), int(Y)
            assert X * X - D * Y * Y == clearing * clearing * sign**r
            assert gcd(X + a * Y, b * Y) <= abs(b) * clearing
            norm = (X + a * Y) ** 2 + (b * Y) ** 2
            assert norm == clearing * clearing * (xy[2 * r][0] + a * xy[2 * r][1])
            norms.append(norm)
            if D == 5 and (a, b, sign) == (1, -2, -1):
                assert (X + a * Y, b * Y) == (
                    2 * fibonacci(r + 1), -2 * fibonacci(r)
                )

        delta = clearing**4 * b * b * xy[2][1] ** 2
        assert delta.denominator == 1 and delta > 0
        delta = int(delta)
        for r in range(len(norms) - 2):
            assert norms[r + 2] == recurrence_trace * norms[r + 1] - norms[r]
            assert norms[r] * norms[r + 2] - norms[r + 1] ** 2 == delta
        U_prev, U_curr = 0, 1
        for d in range(1, 7):
            if d > 1:
                U_prev, U_curr = U_curr, recurrence_trace * U_curr - U_prev
            for r in range(len(norms) - d):
                common = gcd(norms[r], norms[r + d])
                assert common * common <= U_curr * U_curr * delta
        fixtures += 1
    print(f"PASS: {fixtures} general norm-+/-1 unit fixtures, clearing and Wronskians")


def check_minus_one_cancellation() -> None:
    for trace in range(1, 8):
        A = trace * trace + 2
        B = A**3 - 3 * A
        first, third = [1, A, 1], [1, B, 1]
        assert weight(convolve(first, third)) == (A + 2) * (B + 2)
        assert B >= 18 and (A + 2) * (B + 2) > 20
    assert convolve([1, 3, 1], [1, -1]) == [1, 2, -2, -1]

    widths = (1, 3, 4, 3, 1)
    rows = (
        (0, 0, 2, 3, 1),
        (0, 1, 4, 1, 0),
        (1, 2, 0, 2, 1),
        (1, 3, 2, 0, 0),
    )
    assert sum(widths) == 12
    assert tuple(max(row[s] for row in rows) - min(row[s] for row in rows) for s in range(5)) == widths
    assert {sum(row) for row in rows} == {6}

    def powers_as_pairs(middle: int) -> list[tuple[int, int]]:
        powers = [(1, 0), (0, 1)]
        for _ in range(3):
            previous, current = powers[-2], powers[-1]
            powers.append((-previous[0] - middle * current[0],
                           -previous[1] - middle * current[1]))
        return powers

    q_powers = powers_as_pairs(3)
    cubic_powers = powers_as_pairs(18)
    def value(row: tuple[int, ...], powers: list[tuple[int, int]]) -> tuple[int, int]:
        return (sum(x * coefficient for (x, _), coefficient in zip(powers, row)),
                sum(y * coefficient for (_, y), coefficient in zip(powers, row)))

    assert {value(row, q_powers) for row in rows} == {(-1, -3)}
    assert len({value(row, cubic_powers) for row in rows}) == len(rows)

    def multiply(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
        a, b = left
        c, d = right
        return a * c - b * d, a * d + b * c

    for n in range(4, 13):
        points = []
        for row in rows:
            point = (1, 0)
            for s, exponent in enumerate(row):
                block = (fibonacci(n + s + 1), fibonacci(n + s))
                for _ in range(exponent):
                    point = multiply(point, block)
                for _ in range(widths[s] - exponent):
                    point = multiply(point, (block[0], -block[1]))
            points.append(point)
        assert len(set(points)) == 4
        assert len({x * x + y * y for x, y in points}) == 1
    print("PASS: norm-minus-one cancellation polynomials and four-row K=12 fixture")


if __name__ == "__main__":
    check_fibonacci_orbit()
    check_low_height_multiples()
    check_general_quadratic_units()
    check_minus_one_cancellation()
    print("PASS: finer Pell recurrence, norms, shifted gcds, parity, and resonance")
