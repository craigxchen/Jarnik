"""Exact checks of the cotangent spectrum and extreme projectors."""

from fractions import Fraction
from math import factorial, gcd, isqrt, lcm
from itertools import combinations, combinations_with_replacement

from check_integer_cotangent_normalization import edge_quotient, check_clique, primitive_tuple
from check_least_radius_formula import gcd_all, exact_div, lcm_all


ZERO = (0, 0)
ONE = (1, 0)


def add(z, w):
    return z[0] + w[0], z[1] + w[1]


def neg(z):
    return -z[0], -z[1]


def mul(z, w):
    return z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0]


def conj(z):
    return z[0], -z[1]


def div(z, w):
    denominator = w[0]*w[0] + w[1]*w[1]
    numerator = mul(z, conj(w))
    return (Fraction(numerator[0], denominator),
            Fraction(numerator[1], denominator))


def power(z, exponent):
    ans = ONE
    for _ in range(exponent):
        ans = mul(ans, z)
    return ans


def identity(n):
    return [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]


def matrix_mul(a, b):
    n = len(a)
    ans = [[ZERO for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                ans[i][j] = add(ans[i][j], mul(a[i][k], b[k][j]))
    return ans


def shifted(b, eigenvalue):
    ans = [row[:] for row in b]
    for i in range(len(b)):
        ans[i][i] = add(ans[i][i], neg(eigenvalue))
    return ans


def primitive_vector(vector):
    denominator = lcm(*(Fraction(c).denominator for z in vector for c in z))
    cleared = [tuple(int(c*denominator) for c in z) for z in vector]
    common = gcd_all(cleared)
    primitive = [exact_div(z, common) for z in cleared]
    multiplier = div((denominator, 0), common)
    return primitive, multiplier


def dot(u, v):
    result = ZERO
    for z, w in zip(u, v):
        result = add(result, mul(z, w))
    return result


def supported_on_scale(integer, scale):
    assert integer > 0
    while (common := gcd(integer, 2*scale)) > 1:
        integer //= common
    return integer == 1


def check_spectral(xs, scale):
    # Also retain the actual least radius, rather than testing a relaxed matrix.
    check_clique(xs, scale)
    n = len(xs) + 1
    q = [[0]*n for _ in range(n)]
    for i, x in enumerate(xs, 1):
        q[0][i], q[i][0] = x, -x
    for i, x in enumerate(xs, 1):
        for j, y in enumerate(xs, 1):
            if i < j:
                q[i][j] = edge_quotient(x, y, scale)
                q[j][i] = -q[i][j]
    row_sums = [sum(row) for row in q]
    b = [[(row_sums[i], 0) if i == j else (q[i][j], scale)
          for j in range(n)] for i in range(n)]
    nodes = [ONE] + [div((x, -scale), (x, scale)) for x in xs]
    weights = []
    for i in range(n):
        product = ONE
        for j in range(n):
            if i != j:
                product = mul(product, add(nodes[i], neg(nodes[j])))
        weights.append(div(ONE, product))

    eigenvalues = [(0, scale*(2*r-(n-1))) for r in range(n)]
    for r, eigenvalue in enumerate(eigenvalues):
        vector = [mul(weights[i], power(nodes[i], r)) for i in range(n)]
        for i in range(n):
            value = ZERO
            for j in range(n):
                value = add(value, mul(b[i][j], vector[j]))
            assert value == mul(eigenvalue, vector[i])

    # The distinct-node Vandermonde makes these n eigenvectors independent.
    assert len(set(nodes)) == n

    for i in range(n):
        for j in range(n):
            assert b[i][j] == div(mul(nodes[i], conj(b[i][j])), nodes[j])
            hermitian = add(b[i][j], conj(b[j][i]))
            assert hermitian == ((2*row_sums[i], 0) if i == j else ZERO)

    frobenius_squared = sum(z[0]**2+z[1]**2 for row in b for z in row)
    eigen_norm_squared = sum(z[0]**2+z[1]**2 for z in eigenvalues)
    assert frobenius_squared-eigen_norm_squared == 2*sum(s*s for s in row_sums)

    product_nodes = ONE
    for node in nodes:
        product_nodes = mul(product_nodes, node)
    common_denominator = mul(power((0, 2*scale), n-1), (factorial(n-1), 0))
    for r in (0, n-1):
        numerator = identity(n)
        denominator = ONE
        for s, eigenvalue in enumerate(eigenvalues):
            if s != r:
                numerator = matrix_mul(numerator, shifted(b, eigenvalue))
                denominator = mul(denominator,
                                  add(eigenvalues[r], neg(eigenvalue)))
        for i in range(n):
            for j in range(n):
                projector_entry = div(numerator[i][j], denominator)
                if r == n-1:
                    expected = mul(weights[i], power(nodes[i], n-1))
                else:
                    expected = div(mul(((-1)**(n-1), 0),
                                       mul(product_nodes, weights[i])), nodes[j])
                assert projector_entry == expected
                cleared = mul(common_denominator, projector_entry)
                assert all(Fraction(c).denominator == 1 for c in cleared)

    # The inverse weighted Vandermonde has column j equal to the
    # coefficients of product_(l!=j)(z-r_l).
    inverse_columns = []
    for j in range(n):
        coefficients = [ONE]
        for l in range(n):
            if l == j:
                continue
            updated = [ZERO]*(len(coefficients)+1)
            for k, coefficient in enumerate(coefficients):
                updated[k] = add(updated[k], mul(neg(nodes[l]), coefficient))
                updated[k+1] = add(updated[k+1], coefficient)
            coefficients = updated
        inverse_columns.append(coefficients)

    right_vectors, left_vectors, multipliers, pairings = [], [], [], []
    for r in range(n):
        raw = [mul(weights[i], power(nodes[i], r)) for i in range(n)]
        right, multiplier = primitive_vector(raw)
        left, _ = primitive_vector([inverse_columns[j][r] for j in range(n)])
        right_vectors.append(right)
        left_vectors.append(left)
        multipliers.append(multiplier)
        pairing = dot(left, right)
        pairings.append(pairing)
        assert supported_on_scale(pairing[0]**2+pairing[1]**2, scale)
        gap_product = ONE
        for s in range(n):
            if s != r:
                gap_product = mul(gap_product,
                                  add(eigenvalues[r], neg(eigenvalues[s])))
        quotient = div(gap_product, pairing)
        assert all(Fraction(c).denominator == 1 for c in quotient)

    for r in range(n):
        for s in range(n):
            if r != s:
                assert dot(left_vectors[r], right_vectors[s]) == ZERO

    # Exact determinant from the Vandermonde formula, followed by the
    # actual primitive Gaussian multiplier of every column.
    determinant = ONE
    for i in range(n):
        determinant = mul(determinant, weights[i])
        determinant = mul(determinant, multipliers[i])
        for j in range(i+1, n):
            determinant = mul(determinant, add(nodes[j], neg(nodes[i])))
    assert all(Fraction(c).denominator == 1 for c in determinant)
    determinant_norm = int(determinant[0]**2+determinant[1]**2)
    assert determinant_norm > 0
    assert supported_on_scale(determinant_norm, scale), (xs, scale, determinant)
    pairing_product = ONE
    for pairing in pairings:
        pairing_product = mul(pairing_product, pairing)
    left_determinant = div(pairing_product, determinant)
    assert all(Fraction(c).denominator == 1 for c in left_determinant)
    assert supported_on_scale(int(left_determinant[0]**2+left_determinant[1]**2), scale)

    # Exact individual column heights, with conductor integers recovered
    # from subset products rather than from rational-prime factorization.
    circle_rows = primitive_tuple(xs, scale)
    spread_integers = {0: 1}
    for m in range(1, n//2+1):
        products = []
        for indices in combinations(range(n), m):
            product = ONE
            for i in indices:
                product = mul(product, circle_rows[i])
            products.append(product)
        ratio = exact_div(lcm_all(products), gcd_all(products))
        ratio_norm = ratio[0]**2+ratio[1]**2
        spread = isqrt(ratio_norm)
        assert spread*spread == ratio_norm
        spread_integers[m] = spread
    raw_norm_squared = sum(w[0]**2+w[1]**2 for w in weights)
    for r, multiplier in enumerate(multipliers):
        m = n-r-1
        spread = spread_integers[min(m, n-m)]
        defect = Fraction(spread*(multiplier[0]**2+multiplier[1]**2))
        assert defect.denominator == 1
        assert (2*scale)**(2*(n-1)) % defect.numerator == 0
        column_norm_squared = sum(z[0]**2+z[1]**2 for z in right_vectors[r])
        assert column_norm_squared == raw_norm_squared*defect/spread
    if n % 2 == 0:
        h = n//2
        central = right_vectors[h-1]
        assert all(z[1] == 0 for z in central) or all(z[0] == 0 for z in central)
        coordinate = 0 if all(z[1] == 0 for z in central) else 1
        coefficients = [z[coordinate] for z in central]
        for degree in range(h):
            moment = ZERO
            for c, z in zip(coefficients, circle_rows):
                moment = add(moment, mul((c, 0), power(z, degree)))
            assert moment == ZERO
    else:
        h = n//2
        assert any(all(mul(unit, conj(a)) == b
                       for a, b in zip(right_vectors[h-1], right_vectors[h]))
                   for unit in (ONE, (-1, 0), (0, 1), (0, -1)))


def main():
    cases = [((x,), 1) for x in (1, 17, 1234)]
    cases.extend([
        ((3723, 3456, 7184), 24),
        ((282, 316, 342, 542), 6),
        ((-18, -6, -3, 0, 2, 6, 12), 6),
    ])
    for size in range(2, 8):
        scale = lcm(*range(1, size))
        cases.append((tuple(scale*(20+j) for j in range(size)), scale))
    for xs, scale in cases:
        check_spectral(xs, scale)
    profiles = 0
    for n in range(2, 9):
        for exponents in combinations_with_replacement(range(-2, 3), n):
            column_minima = []
            for r in range(n):
                values = [(r-n+i+1)*exponents[i]-sum(exponents[:i])
                          for i in range(n)]
                minimum = -sum(exponents[:n-r-1])
                assert min(values) == minimum
                column_minima.append(minimum)
            assert sum(column_minima) == -sum((n-i-1)*e
                                              for i, e in enumerate(exponents))
            profiles += 1
    print("PASS:", len(cases), "actual cliques; exact spectrum, conjugacy,")
    print("extreme projectors, integral denominator, and nonnormality identity.")
    print("PASS: primitive eigenvector pairings divide spectral gaps;")
    print("primitive eigenbasis determinants have prime support in 2L.")
    print("PASS: exact primitive column heights, bounded integer defects,")
    print("central real moment relations and conjugate-column pairs.")
    print("PASS:", profiles, "ordered valuation profiles, including negative levels.")


if __name__ == "__main__":
    main()
