"""Exact certificates for shrinking_arcs_with_bounded_pair_residues.md."""

from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, prod


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def conj(z):
    return z[0], -z[1]


def norm(z):
    return z[0] ** 2 + z[1] ** 2


def div_exact(z, w):
    a, b = mul(z, conj(w))
    q = norm(w)
    assert a % q == b % q == 0
    return a // q, b // q


def ggcd(z, w):
    while w != (0, 0):
        a, b = mul(z, conj(w))
        q = norm(w)
        nearest = ((2 * a + q) // (2 * q), (2 * b + q) // (2 * q))
        p = mul(nearest, w)
        z, w = w, (z[0] - p[0], z[1] - p[1])
    return z


def glcm(rows):
    value = (1, 0)
    for row in rows:
        value = mul(div_exact(value, ggcd(value, row)), row)
    return value


def permanent(matrix):
    size = len(matrix)
    return sum((prod(matrix[i][s[i]] for i in range(size))
                for s in permutations(range(size))), Fraction(0))


def check_remainder(rows):
    size = len(rows)
    norms = [norm(z) for z in rows]
    dots = [[mul(z, conj(w))[0] for w in rows] for z in rows]
    assert all(value > 0 for row in dots for value in row)
    per_s = prod(norms) * permanent([[Fraction(1, x) for x in row]
                                   for row in dots])
    linear = sum((Fraction(mul(rows[i], conj(rows[j]))[1] ** 2, norms[i] * norms[j])
                  for i, j in combinations(range(size), 2)), Fraction(0))
    remainder = per_s - factorial(size) - factorial(size - 1) * linear
    assert remainder > 0
    assert remainder.denominator * remainder >= 1
    return remainder.denominator.bit_length()


def primary_unit(z):
    for unit in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        candidate = mul(z, conj(unit))
        difference = (candidate[0] - 1, candidate[1])
        numerator = mul(difference, (-2, -2))
        if numerator[0] % 8 == numerator[1] % 8 == 0:
            return unit
    raise AssertionError("odd-norm Gaussian integer has no primary associate")


def main():
    r, b = 1, 0
    denominator_bits = []
    cases = 0
    for index in range(1, 7):
        r, b = 3 * r + 4 * b, 2 * r + 3 * b
        assert r * r - 2 * b * b == 1
        a = r + b
        g = (a, -b)
        d = norm(g)
        n = a * a - b * b + 2 * a * b
        assert mul((1, 1), mul(g, g)) == (n, 1)
        assert n * n + 1 == 2 * d * d
        direction = mul((1, -1), conj(g))
        all_rows = [(g[0] + j * direction[0], g[1] + j * direction[1])
                    for j in range(20)]
        for size in range(2, 11):
            rows = all_rows[:size]
            for j, row in enumerate(rows):
                assert mul((1, 1), mul(g, row)) == (n + 2 * j * d, 1)
                assert norm(row) % 2 == 1
                assert norm(ggcd(row, conj(row))) == 1
            residue_product = 1
            for i, j in combinations(range(size), 2):
                assert mul(rows[i], conj(rows[j]))[1] == j - i
                common_norm = norm(ggcd(rows[i], rows[j]))
                assert (j - i) % common_norm == 0
                residue_product *= (j - i) // common_norm
            differences = prod(j - i for i, j in combinations(range(size), 2))
            assert 1 <= residue_product <= differences
            lcd = glcm(rows)
            points = [div_exact(mul(conj(lcd), row), conj(row)) for row in rows]
            assert len(set(points)) == size
            assert all(norm(z) == norm(lcd) for z in points)
            common = points[0]
            for z in points[1:]:
                common = ggcd(common, z)
            assert norm(common) == 1
            product_norms = prod(norm(z) for z in rows)
            assert product_norms <= norm(lcd) * differences
            assert norm(lcd) <= product_norms
            cases += 1
        denominator_bits.append(check_remainder(all_rows[:3]))
        check_remainder(all_rows[:4])
        even_rows = all_rows[::2]
        assert all(primary_unit(z) == primary_unit(g) for z in even_rows)
        assert norm(ggcd(even_rows[0], even_rows[1])) == 1
        for i, j in combinations(range(len(even_rows)), 2):
            assert mul(even_rows[i], conj(even_rows[j]))[1] == 2 * (j - i)
        check_remainder(even_rows[:3])
    print(f"PASS: {cases} primitive Gaussian circles with 2 through 10 rows.")
    print("Exact pencil determinants, bounded pair residues, and full-radius bounds pass.")
    print("Three-point first-remainder reduced-denominator bit lengths:", denominator_bits)
    print("Three- and four-point positive rational Taylor remainders pass at six Pell indices.")
    print("Even-index selection has one fixed primary unit and the same remainder property.")
    print("The m>4 families are not endpoint counterexamples.")


if __name__ == "__main__":
    main()
