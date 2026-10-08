"""Finite checks for the first-smaller-cut Chow basepoint note.

This checks matching-space coordinates, small-field base loci, sample
birational ratios, small-n degree counts, and the exact m=400 inequality.
It does not compute a Chow form or verify the Gaussian arithmetic theorem.
"""

from fractions import Fraction
from itertools import combinations, product
from math import comb


def matching_value(z, a, b, c, d):
    return (z[a] - z[b]) * (z[c] - z[d])


def matchings(n):
    for a, b, c, d in combinations(range(n), 4):
        yield a, b, c, d
        yield a, c, b, d
        yield a, d, b, c


def basis_pairs(n):
    return [(i, j) for i, j in combinations(range(n - 1), 2)
            if (i, j) != (0, 1)]


def basis_value(z, i, j):
    # The quotient-by-translation chart is z_(n-1)=0.
    return z[i] * z[j] - z[0] * z[1]


def rank_mod_p(rows, prime=1_000_003):
    pivots = {}
    for row in rows:
        value = [entry % prime for entry in row]
        for col in sorted(pivots):
            if value[col]:
                scalar = value[col]
                value = [(a - scalar * b) % prime
                         for a, b in zip(value, pivots[col])]
        col = next((j for j, entry in enumerate(value) if entry), None)
        if col is not None:
            inverse = pow(value[col], -1, prime)
            pivots[col] = [(entry * inverse) % prime for entry in value]
    return len(pivots)


def check_matching_space():
    for n in range(5, 10):
        pairs = list(combinations(range(n - 1), 2))
        rows = []
        for i, j in basis_pairs(n):
            rows.append([int((a, b) == (i, j)) - int((a, b) == (0, 1))
                         for a, b in pairs])
        t = n * (n - 3) // 2
        assert len(rows) == t
        assert rank_mod_p(rows) == t
        assert comb(n - 1, 2) - 1 == t

        # The n projective frame points in this chart are the coordinate
        # vectors and the all-minus-one vector.
        frames = []
        for i in range(n - 1):
            z = [0] * n
            z[i] = 1
            frames.append(z)
        frames.append([-1] * (n - 1) + [0])
        for z in frames:
            assert all(basis_value(z, i, j) == 0
                       for i, j in basis_pairs(n))


def is_frame_point(z):
    # The nonzero quotient vector has one value on n-1 labels and a
    # different value on the remaining label.
    return any(len({z[j] for j in range(len(z)) if j != i}) == 1
               and z[i] != z[(i + 1) % len(z)]
               for i in range(len(z)))


def check_small_field_base_locus():
    for n in range(5, 9):
        matching_indices = list(matchings(n))
        count = 0
        for prefix in product(range(3), repeat=n - 1):
            z = prefix + (0,)
            if all(value == 0 for value in z):
                continue
            is_base = all(matching_value(z, *edges) % 3 == 0
                          for edges in matching_indices)
            assert is_base == is_frame_point(z), (n, z)
            count += is_base
        # A projective frame point has two representatives in F_3.
        assert count == 2 * n, (n, count)


def check_birational_ratios():
    for n in range(5, 10):
        z = tuple(3 * j * j + 2 * j + 1 for j in range(n))
        for e in range(2, n):
            a, b = next(combinations((j for j in range(n)
                                      if j not in (0, 1, e)), 2))
            numerator = matching_value(z, a, b, 0, 1)
            denominator = matching_value(z, a, b, 0, e)
            ratio = Fraction(numerator, denominator)
            recovered = z[0] - Fraction(z[0] - z[1], ratio)
            assert recovered == z[e]


def check_degree_examples():
    expected = {5: 3, 6: 10, 7: 25, 8: 56, 9: 119, 10: 246}
    for n, degree in expected.items():
        dimension = n - 2
        # On the blowup of P^dimension at n points, the only nonzero
        # top products are H^dimension=1 and E_i^dimension=(-1)^(dimension-1).
        intersection = 2 ** dimension + n * (-1) ** dimension * (-1) ** (dimension - 1)
        assert intersection == degree == 2 ** (n - 2) - n


def coefficient_of_trinomial_power(m):
    coefficients = [1]
    for _ in range(m):
        next_coefficients = [0] * (len(coefficients) + 2)
        for j, value in enumerate(coefficients):
            next_coefficients[j] += value
            next_coefficients[j + 1] += value
            next_coefficients[j + 2] += value
        coefficients = next_coefficients
    return coefficients


def check_m400_inequality():
    m = 400
    q = m // 2
    n = q + 1
    coefficients = coefficient_of_trinomial_power(m)
    r = coefficients[m] - coefficients[m - 1]
    A = m * 2 ** (m - 2) - q * comb(m, q)
    L = (r - 1) // 2
    delta = 2 ** (n - 2) - n
    D = delta * (n - 1)
    assert r > L > 0 and A > 0 and D > 0
    assert 40 * D * A < r - L
    return r, A, L, D


def main():
    check_matching_space()
    check_small_field_base_locus()
    check_birational_ratios()
    check_degree_examples()
    r, A, L, D = check_m400_inequality()
    print("PASS: matching dimensions and frame basepoints for n=5..9;")
    print("      full F_3 matching base loci for n=5..8;")
    print("      birational ratios and blowup degree examples;")
    print("      exact 40*delta*(n-1)*A_400 < r_400-L_400")
    print(f"      (r has {len(str(r))} digits, A has {len(str(A))} digits, "
          f"L has {len(str(L))} digits, D has {len(str(D))} digits).")


if __name__ == "__main__":
    main()
