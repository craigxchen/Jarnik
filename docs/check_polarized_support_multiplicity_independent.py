#!/usr/bin/env python3
"""Bounded exhaustive audit of polarization, grid distance, and cut width."""

from itertools import product
from math import comb


def exponent_box(bounds):
    return list(product(*(range(bound + 1) for bound in bounds)))


def grid_points(bounds):
    return exponent_box(bounds)


def polynomial_value(coefficients, exponents, point, prime):
    return sum(
        coefficient
        * product_value(pow(x, power, prime) for x, power in zip(point, exponent))
        for coefficient, exponent in zip(coefficients, exponents)
    ) % prime


def product_value(values):
    answer = 1
    for value in values:
        answer *= value
    return answer


def translated_coefficients(coefficients, exponents, prime):
    """Coefficients after substituting z_i=1+y_i."""
    bounds = tuple(max(exponent[i] for exponent in exponents) for i in range(len(exponents[0])))
    translated = {}
    for gamma in exponent_box(bounds):
        value = 0
        for coefficient, exponent in zip(coefficients, exponents):
            if all(g <= a for g, a in zip(gamma, exponent)):
                value += coefficient * product_value(
                    comb(a, g) for a, g in zip(exponent, gamma)
                )
        translated[gamma] = value % prime
    return translated


def multiplicity_at_one(coefficients, exponents, prime):
    translated = translated_coefficients(coefficients, exponents, prime)
    nonzero_degrees = [sum(exponent) for exponent, value in translated.items() if value]
    return min(nonzero_degrees) if nonzero_degrees else None


def subsets(size):
    return range(1 << size)


def popcount(mask):
    return bin(mask).count("1")


def group_count(mask, group):
    return sum((mask >> position) & 1 for position in group)


def check_polarization_and_indicator():
    """Exhaust all 3^9 bidegree-(2,2) polynomials over F_3."""
    prime = 3
    bounds = (2, 2)
    exponents = exponent_box(bounds)
    groups = ((0, 1), (2, 3))
    variable_count = 4
    checked = 0
    for coefficients in product(range(prime), repeat=len(exponents)):
        if not any(coefficients):
            continue
        coefficient_map = dict(zip(exponents, coefficients))
        multiplicity = multiplicity_at_one(coefficients, exponents, prime)

        polarized = {}
        for mask in subsets(variable_count):
            counts = tuple(group_count(mask, group) for group in groups)
            denominator = product_value(
                comb(bound, count) for bound, count in zip(bounds, counts)
            )
            polarized[mask] = coefficient_map[counts] * pow(denominator, -1, prime) % prime

        for exponent in exponents:
            assert sum(
                value
                for mask, value in polarized.items()
                if tuple(group_count(mask, group) for group in groups) == exponent
            ) % prime == coefficient_map[exponent]

        translated_polarized = {}
        for mask in subsets(variable_count):
            translated_polarized[mask] = sum(
                value for larger, value in polarized.items() if mask & ~larger == 0
            ) % prime
        polarized_degrees = [
            popcount(mask) for mask, value in translated_polarized.items() if value
        ]
        assert min(polarized_degrees) == multiplicity

        indicator_coefficients = {}
        for mask in subsets(variable_count):
            indicator_coefficients[mask] = (
                (-1) ** popcount(mask)
                * sum(value for smaller, value in polarized.items() if smaller & ~mask == 0)
            ) % prime
        indicator_degree = max(
            (popcount(mask) for mask, value in indicator_coefficients.items() if value),
            default=-1,
        )
        assert indicator_degree <= variable_count - multiplicity
        for mask in subsets(variable_count):
            for group in groups:
                swapped = mask
                left, right = group
                left_bit, right_bit = (mask >> left) & 1, (mask >> right) & 1
                if left_bit != right_bit:
                    swapped ^= (1 << left) | (1 << right)
                assert indicator_coefficients[mask] == indicator_coefficients[swapped]

        for vertex in subsets(variable_count):
            value = sum(
                coefficient
                for mask, coefficient in indicator_coefficients.items()
                if mask & ~vertex == 0
            ) % prime
            expected = (-1) ** popcount(vertex) * polarized[vertex] % prime
            assert value == expected
            counts = tuple(group_count(vertex, group) for group in groups)
            assert (value != 0) == (coefficient_map[counts] != 0)
        checked += 1
    assert checked == 3**9 - 1


def check_grid_instance(bounds, prime):
    exponents = exponent_box(bounds)
    points = grid_points(bounds)
    corners = list(product(*((0, bound) for bound in bounds)))
    checked = 0
    for coefficients in product(range(prime), repeat=len(exponents)):
        if not any(coefficients):
            continue
        support = [
            point
            for point in points
            if polynomial_value(coefficients, exponents, point, prime)
        ]
        assert support
        total_degree = max(
            sum(exponent)
            for coefficient, exponent in zip(coefficients, exponents)
            if coefficient
        )
        distance_sum = sum(
            min(sum(abs(x - y) for x, y in zip(corner, point)) for point in support)
            for corner in corners
        )
        assert 2 * distance_sum <= total_degree * len(corners)
        checked += 1
    return checked


def check_grid_lemma():
    # These exhaust arbitrary polynomials, not only polarized examples.
    assert check_grid_instance((2, 2), 3) == 3**9 - 1
    assert check_grid_instance((1, 1, 1), 3) == 3**8 - 1


def cut_representatives(variable_count):
    answer = []
    for mask in range(1, (1 << variable_count) - 1):
        size = popcount(mask)
        if size < variable_count - size:
            answer.append(mask)
        elif size == variable_count - size and mask & 1:
            answer.append(mask)
    return answer


def homogeneous_exponents(bounds):
    degree = sum(bounds) // 2
    return [exponent for exponent in exponent_box(bounds) if sum(exponent) == degree]


def reciprocal_index(exponent, bounds):
    return tuple(bound - power for bound, power in zip(bounds, exponent))


def reciprocal_coefficient_vectors(bounds, prime):
    exponents = homogeneous_exponents(bounds)
    index = {exponent: position for position, exponent in enumerate(exponents)}
    orbits = []
    seen = set()
    for exponent in exponents:
        if exponent in seen:
            continue
        partner = reciprocal_index(exponent, bounds)
        seen.add(exponent)
        seen.add(partner)
        orbits.append((index[exponent], index[partner]))

    choices = []
    for left, right in orbits:
        if left == right:
            choices.append([(value,) for value in range(prime)])
        else:
            choices.append([(0, 0)] + list(product(range(1, prime), repeat=2)))
    for orbit_values in product(*choices):
        coefficients = [0] * len(exponents)
        for (left, right), values in zip(orbits, orbit_values):
            if left == right:
                coefficients[left] = values[0]
            else:
                coefficients[left], coefficients[right] = values
        if any(coefficients):
            yield exponents, tuple(coefficients)


def doubled_width(exponents, coefficients, bounds):
    support = [exponent for exponent, coefficient in zip(exponents, coefficients) if coefficient]
    answer = 0
    degree = sum(bounds) // 2
    for cut in cut_representatives(len(bounds)):
        height = max(
            sum(2 * exponent[i] - bounds[i] for i in range(len(bounds)) if cut >> i & 1)
            for exponent in support
        )
        corner = tuple(bounds[i] if cut >> i & 1 else 0 for i in range(len(bounds)))
        distance = min(
            sum(abs(x - y) for x, y in zip(corner, exponent)) for exponent in support
        )
        assert distance == degree - height
        answer += height
    return answer


def check_reciprocal_width(bounds, prime, expected_count):
    checked = 0
    factor = 1 << (len(bounds) - 2)  # doubled form of 2^(n-3)*m
    for exponents, coefficients in reciprocal_coefficient_vectors(bounds, prime):
        multiplicity = multiplicity_at_one(coefficients, exponents, prime)
        assert multiplicity is not None
        assert doubled_width(exponents, coefficients, bounds) >= factor * multiplicity
        checked += 1
    assert checked == expected_count


def check_central_symmetry_normalization():
    # Arbitrary paired nonzero coefficients are allowed; only support symmetry
    # is imposed.  The bounds are small enough for complete enumeration.
    check_reciprocal_width((1, 1, 1, 1), 5, 17**3 - 1)
    check_reciprocal_width((2, 2, 2), 5, 5 * 17**3 - 1)


def main():
    check_polarization_and_indicator()
    check_grid_lemma()
    check_central_symmetry_normalization()
    print(
        "PASS: exhaustive small polarization, grid-distance, and reciprocal "
        "cut-width checks found no counterexample."
    )


if __name__ == "__main__":
    main()
