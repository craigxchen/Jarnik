"""Exact certificates for the normalized Vieta bound and moving frame.

Only the standard library is needed.  The literal full-cut fixture checks
divisibility and content, not endpoint angles or negligible residues.  The
moving-frame example likewise lies outside the full-profile hypotheses.
"""

from fractions import Fraction
from itertools import combinations, permutations, product
from math import gcd, isqrt, prod


def poly_add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for exponent, coefficient in polynomial.items():
            result[exponent] = result.get(exponent, 0) + coefficient
    return {exponent: coefficient for exponent, coefficient in result.items()
            if coefficient}


def poly_scale(polynomial, factor):
    return {exponent: factor * coefficient
            for exponent, coefficient in polynomial.items() if factor * coefficient}


def poly_mul(first, second):
    result = {}
    for e, a in first.items():
        for f, b in second.items():
            exponent = tuple(x + y for x, y in zip(e, f))
            result[exponent] = result.get(exponent, 0) + a * b
    return {exponent: coefficient for exponent, coefficient in result.items()
            if coefficient}


def square(polynomial):
    return poly_mul(polynomial, polynomial)


def check_coefficient_norm_identity():
    variables = []
    for i in range(4):
        exponent = tuple(int(i == j) for j in range(4))
        variables.append({exponent: 1})
    a, b, c, d = variables
    ac, ad = poly_mul(a, c), poly_mul(a, d)
    bc, bd = poly_mul(b, c), poly_mul(b, d)
    middle = poly_add(ad, bc)
    coefficient_norm = poly_add(square(ac), square(middle), square(bd))
    factor_norm_product = poly_mul(poly_add(square(a), square(b)),
                                   poly_add(square(c), square(d)))
    residual = poly_add(poly_scale(coefficient_norm, 2),
                        poly_scale(factor_norm_product, -1),
                        poly_scale(square(poly_add(ac, bd)), -1),
                        poly_scale(square(middle), -1))
    assert residual == {}

    vectors = [(a, b) for a, b in product(range(-6, 7), repeat=2)
               if gcd(a, b) == 1]
    cases = 0
    infinity_cases = 0
    for (a, b), (c, d) in product(vectors, repeat=2):
        coefficients = (a * c, a * d + b * c, b * d)
        assert gcd(gcd(coefficients[0], coefficients[1]), coefficients[2]) == 1
        coefficient_norm = sum(x * x for x in coefficients)
        factor_norm_product = (a * a + b * b) * (c * c + d * d)
        assert (2 * coefficient_norm - factor_norm_product
                == (a * c + b * d) ** 2 + (a * d + b * c) ** 2)
        # Squared, exponentiated form of the logarithmic height inequality.
        assert (6 * max(abs(x) for x in coefficients) ** 2
                >= max(abs(a), abs(b)) ** 2 * max(abs(c), abs(d)) ** 2)
        infinity_cases += int(a == 0 or c == 0)
        cases += 1
    assert cases == 9216 and infinity_cases > 0
    return cases


def det(first, second):
    return first[0] * second[1] - first[1] * second[0]


def matrix_vector(matrix, vector):
    return tuple(sum(a * b for a, b in zip(row, vector)) for row in matrix)


def norm(vector):
    return sum(x * x for x in vector)


def max_height(vector):
    return max(abs(x) for x in vector)


def primitive(vector):
    content = gcd(*vector)
    assert content > 0
    return tuple(x // content for x in vector), content


def internal(cut, u, v):
    return int(bool(cut & (1 << u)) and bool(cut & (1 << v)))


def matching_exponents(cut, a, b, c, i):
    return (internal(cut, c, a) + internal(cut, i, b),
            internal(cut, b, c) + internal(cut, a, i),
            internal(cut, a, b) + internal(cut, c, i))


def check_cut_counts():
    cases = 0
    for a, b, c, i in permutations(range(5), 4):
        d_weight = 0
        matching_weight = 0
        matching_totals = [0, 0, 0]
        for cut in range(32):
            d_exponent = min(internal(cut, b, c), internal(cut, c, a))
            assert d_exponent == int(all(cut & (1 << j) for j in (a, b, c)))
            exponents = matching_exponents(cut, a, b, c, i)
            common = min(exponents)
            # Any two matching products have the full three-product minimum.
            assert all(min(exponents[j], exponents[k]) == common
                       for j, k in combinations(range(3), 2))
            inside = sum(bool(cut & (1 << j)) for j in (a, b, c, i))
            assert common == max(0, inside - 2)
            assert min(exponents[0] - common, exponents[1] - common) == 0
            d_weight += d_exponent
            matching_weight += common
            matching_totals = [x + y for x, y in zip(matching_totals, exponents)]
            cases += 1
        assert d_weight == 4
        assert matching_weight == 12
        assert matching_totals == [16, 16, 16]
        row_weight = sum(bool(cut & (1 << a)) for cut in range(32))
        bracket_weight = sum(internal(cut, a, b) for cut in range(32))
        assert row_weight == 16 and bracket_weight == 8
        assert bracket_weight + row_weight // 2 - d_weight == 12  # h(M)
        assert 3 * bracket_weight - 2 * d_weight == 16  # det(M)
        assert matching_totals[0] - matching_weight == 4  # normalized row
        assert 3 * bracket_weight - d_weight - matching_weight == 8  # content
    assert cases == 3840
    return cases


def gaussian_mul(first, second):
    a, b = first
    c, d = second
    return a * c - b * d, a * d + b * c


def literal_full_cut_fixture():
    blocks = []
    used_norms = set()
    for a in range(1, 40):
        for b in range(1, a + 1):
            p = a * a + b * b
            if (p > 2 and p not in used_norms
                    and all(p % d for d in range(2, isqrt(p) + 1))):
                blocks.append((a, b))
                used_norms.add(p)
            if len(blocks) == 31:
                break
        if len(blocks) == 31:
            break
    assert len(blocks) == 31
    rows = []
    for i in range(5):
        value = (1, 0)
        for cut, block in enumerate(blocks, 1):
            if cut & (1 << i):
                value = gaussian_mul(value, block)
        assert gcd(*value) == 1 and norm(value) % 2 == 1
        rows.append(value)
    return rows, {cut: norm(block) for cut, block in enumerate(blocks, 1)}


def check_frame_fixture(rows, cut_norms=None):
    assert len(rows) == 5 and all(gcd(*row) == 1 for row in rows)
    brackets = {(u, v): det(rows[u], rows[v])
                for u in range(5) for v in range(5)}
    assert all(brackets[u, v] for u, v in combinations(range(5), 2))
    corrections = {}
    if cut_norms is not None:
        for u, v in combinations(range(5), 2):
            core = prod(n for cut, n in cut_norms.items() if internal(cut, u, v))
            assert abs(brackets[u, v]) % core == 0
            corrections[u, v] = corrections[v, u] = abs(brackets[u, v]) // core

    cases = 0
    for a, b, c, i in permutations(range(5), 4):
        A, B, C = brackets[b, c], brackets[c, a], brackets[a, b]
        d = gcd(A, B)
        matrix = tuple((A * rows[a][j] // d, B * rows[b][j] // d) for j in range(2))
        assert gcd(gcd(*matrix[0]), gcd(*matrix[1])) == 1
        determinant = det(matrix[0], matrix[1])
        assert A * B * C % (d * d) == 0
        assert determinant == A * B * C // (d * d)
        assert C % d == 0
        assert matrix_vector(matrix, (1, 1)) == tuple(-C // d * x for x in rows[c])

        raw = (B * brackets[i, b], A * brackets[a, i])
        v, s = primitive(raw)
        third = C * brackets[c, i]
        assert gcd(s, third) == s
        assert gcd(raw[0], third) == gcd(raw[1], third) == s
        image = matrix_vector(matrix, v)
        assert A * B * C % (d * s) == 0
        scale = A * B * C // (d * s)
        assert image == tuple(scale * x for x in rows[i])
        reduced, content = primitive(image)
        assert content == abs(scale)
        assert abs(determinant) % content == 0
        assert norm(image) == content * content * norm(rows[i])
        stretch_squared = Fraction(norm(image), norm(v))
        assert Fraction(norm(v)) * stretch_squared / (content * content) == norm(rows[i])

        matrix_height = max(abs(x) for row in matrix for x in row)
        assert max_height(reduced) <= 2 * matrix_height * max_height(v)
        assert max_height(v) <= 2 * matrix_height * max_height(reduced)

        if cut_norms is not None:
            d_core = prod(n for cut, n in cut_norms.items()
                          if all(cut & (1 << j) for j in (a, b, c)))
            s_core = prod(n ** min(matching_exponents(cut, a, b, c, i))
                          for cut, n in cut_norms.items())
            assert d % d_core == 0 and s % s_core == 0
            assert corrections[b, c] * corrections[c, a] % (d // d_core) == 0
            correction_product = prod(corrections[u, v] for u, v in
                                      ((c, a), (i, b), (b, c), (a, i)))
            assert correction_product % (s // s_core) == 0
        cases += 1
    return cases


def check_moving_frame():
    old_row, new_row = (12, 1), (22, 1)
    U = ((0, 1), (-1, 12))
    V = ((1, 0), (2, 1))
    expected_norms = (865, 53065, 50003000065)
    for N, expected in zip((1, 10, 10000), expected_norms):
        S = ((1, N), (0, 1))
        assert all(det(matrix[0], matrix[1]) == 1 for matrix in (U, S, V))
        old = matrix_vector(V, matrix_vector(S, matrix_vector(U, old_row)))
        new = matrix_vector(V, matrix_vector(S, matrix_vector(U, new_row)))
        assert old == (1, 2)
        assert new == (1 - 10 * N, -8 - 20 * N)
        assert gcd(*new) == 1 and norm(new) % 2 == 1
        assert det(old, new) == -10
        assert norm(old) == 5 and norm(new) == expected
    return len(expected_norms)


def main():
    factors = check_coefficient_norm_identity()
    cuts = check_cut_counts()
    fixtures = [([(4 * j, 1) for j in range(1, 6)], None),
                ([(1, 0), (0, 1), (1, 1), (2, 3), (-2, 5)], None),
                literal_full_cut_fixture()]
    frames = sum(check_frame_fixture(rows, norms) for rows, norms in fixtures)
    moving = check_moving_frame()
    print(f"PASS: polynomial coefficient-norm identity and {factors} primitive factor pairs.")
    print(f"PASS: {cuts} labelled cut cases; contents 4 and 12, frame height 12, determinant 16.")
    print(f"PASS: {frames} exact frame/content cases, including a literal 31-block Gaussian fixture.")
    print(f"PASS: {moving} unimodular moving-frame examples (outside the endpoint profile).")


if __name__ == "__main__":
    main()
