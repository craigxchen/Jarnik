"""Exact checks for the local tensor coefficient sharpness construction.

The generic-coalition argument is proved in the accompanying note.
These tests do not construct endpoint circle configurations.
"""

from collections import defaultdict
from itertools import permutations, product
from math import factorial


def sign(perm):
    return (-1) ** sum(perm[i] > perm[j]
                      for i in range(len(perm))
                      for j in range(i + 1, len(perm)))


def determinant_kernel(n, z):
    """Expand det[v | z*v] by actual determinant permutations."""
    coefficients = defaultdict(int)
    for perm in permutations(range(2 * n)):
        value = sign(perm)
        colors = []
        for i, column in enumerate(perm):
            colors.append(column % n)
            if column >= n:
                value *= z[i]
        coefficients[tuple(colors)] += value
    return {key: value for key, value in coefficients.items() if value}


def substituted(poly, z, a, b):
    """Top color derivative divided by the factorial of its degree."""
    output = defaultdict(int)
    for colors, coefficient in poly.items():
        target = list(colors)
        for i, color in enumerate(colors):
            if color == b:
                target[i] = a
                coefficient *= z[i]
        output[tuple(target)] += coefficient
    return {key: value for key, value in output.items() if value}


def determinant_product(n, d, local):
    output = local
    for _ in range(d - 2):
        next_output = {}
        for colors, coefficient in output.items():
            for perm in permutations(range(n)):
                key = colors + perm
                assert key not in next_output
                next_output[key] = coefficient * sign(perm)
        output = next_output
    return output


def bracket_polynomial(first, second):
    """Rows consist of (constant P, linear P, constant Y)."""
    p0, p1, y = first
    q0, q1, v = second
    return (p0 * v - q0 * y, p1 * v - q1 * y)


def coefficient_order(rows, pairs):
    """Multiply bracket polynomials without assuming their order formula."""
    polynomial = [1]
    for i, j in pairs:
        c0, c1 = bracket_polynomial(rows[i], rows[j])
        next_poly = [0] * (len(polynomial) + 1)
        for degree, value in enumerate(polynomial):
            next_poly[degree] += value * c0
            next_poly[degree + 1] += value * c1
        polynomial = next_poly
    return next(i for i, value in enumerate(polynomial) if value)


def matchings(labels):
    if not labels:
        yield ()
        return
    first = labels[0]
    for i in range(1, len(labels)):
        second = labels[i]
        rest = labels[1:i] + labels[i + 1:]
        for tail in matchings(rest):
            yield ((first, second),) + tail


def collision_rows(n, inside):
    # All inside rows tend to infinity; outside rows have distinct slopes.
    return [(0, i + 1, 1) if i in inside else (1, 0, i + 2)
            for i in range(2 * n)]


def main():
    for n in range(2, 5):
        z = [i * i + 2 for i in range(2 * n)]
        local = determinant_kernel(n, z)
        assert len(local) == factorial(2 * n) // (2 ** n)
        for colors, coefficient in local.items():
            expected = 1
            for a in range(n):
                indices = [i for i, color in enumerate(colors) if color == a]
                assert len(indices) == 2
                expected *= z[indices[1]] - z[indices[0]]
            assert abs(coefficient) == expected
        for a in range(n):
            for b in range(n):
                if a != b:
                    assert not substituted(local, z, a, b)
        print(f'n={n}: determinant coefficient and all color-pair checks pass')

    for n, d in ((2, 3), (3, 3), (3, 4)):
        z = [i * i + 2 for i in range(n * d)]
        local = determinant_kernel(n, z[:2 * n])
        polynomial = determinant_product(n, d, local)
        for colors in polynomial:
            assert all(colors.count(a) == d for a in range(n))
        for a in range(n):
            for b in range(n):
                if a != b:
                    assert not substituted(polynomial, z, a, b)
        print(f'(n,d)=({n},{d}): product CB annihilator passes')

    for n in range(2, 6):
        all_pairs = tuple(matchings(tuple(range(2 * n))))
        for s in range(2 * n + 1):
            rows = collision_rows(n, set(range(s)))
            minimum = min(coefficient_order(rows, pairs) for pairs in all_pairs)
            assert minimum == max(0, s - n)

    cases = 0
    for n in range(2, 7):
        chosen_pairs = tuple((2 * a, 2 * a + 1) for a in range(n))
        for d in range(2, 7):
            for occupancies in product(range(d + 1), repeat=n):
                inside = set()
                for a, count in enumerate(occupancies):
                    if count == d:
                        inside.update((2 * a, 2 * a + 1))
                    elif count:
                        inside.add(2 * a)
                rows = collision_rows(n, inside)
                order = coefficient_order(rows, chosen_pairs)
                normalized = order - max(0, len(inside) - n)
                full = occupancies.count(d)
                empty = occupancies.count(0)
                assert normalized == min(full, empty)
                cases += 1
    print(f'PASS: local primitive order and {cases} occupancy witnesses.')


if __name__ == '__main__':
    main()
