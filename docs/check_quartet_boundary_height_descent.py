"""Finite checks for quartet-supported first-smaller height descent.

The symbolic checks use actual rational directions on 10 and 12 rows.
The Gaussian checks use finite independent split-prime blocks on 4 and
6 rows, including correction factors on a core prime. Those rows need
not have a common norm; they test the divisibility mechanism, not a
geometric circle realization. Neither part proves an endpoint height
bound or tests all possible configurations.
"""

from itertools import combinations
from math import gcd, isqrt


def add_terms(*families):
    return [term for family in families for term in family]


def scale_terms(terms, scalar):
    return [(scalar * coefficient, edges) for coefficient, edges in terms]


def multiply_terms(left, right):
    return [(a * b, edges_a + edges_b)
            for a, edges_a in left for b, edges_b in right]


def triangle(a, b, c):
    return [(1, ((a, b), (b, c), (c, a)))]


def quartet_terms():
    X = [(1, ((0, 1), (2, 3)))]
    Z = [(1, ((0, 3), (1, 2)))]
    Q1 = add_terms(scale_terms(multiply_terms(X, X), 3),
                   scale_terms(multiply_terms(X, Z), -1))
    Q2 = add_terms(scale_terms(multiply_terms(X, Z), 3),
                   scale_terms(multiply_terms(Z, Z), -1))
    return Q1, Q2


def graph_value(edges, directions):
    value = 1
    for a, b in edges:
        value *= directions[b] - directions[a]
    return value


def terms_value(terms, directions):
    return sum(coefficient * graph_value(edges, directions)
               for coefficient, edges in terms)


def restrict_graph(edges, inside, m):
    zero = (0,) * m
    polynomial = {zero: 1}
    for a, b in edges:
        if a in inside and b in inside:
            return {}
        if a in inside:
            polynomial = {monomial: -coefficient
                          for monomial, coefficient in polynomial.items()}
        elif b in inside:
            continue
        else:
            updated = {}
            for monomial, coefficient in polynomial.items():
                for label, sign in ((b, 1), (a, -1)):
                    new_monomial = list(monomial)
                    new_monomial[label] += 1
                    key = tuple(new_monomial)
                    updated[key] = updated.get(key, 0) + sign * coefficient
            polynomial = {key: value for key, value in updated.items() if value}
    return polynomial


def restrict_terms(terms, inside, m):
    out = {}
    for scalar, edges in terms:
        for monomial, value in restrict_graph(edges, inside, m).items():
            out[monomial] = out.get(monomial, 0) + scalar * value
    return {monomial: value for monomial, value in out.items() if value}


def quartet_balanced_values(terms):
    values = []
    for i in (1, 2, 3):
        restricted = restrict_terms(terms, {0, i}, 4)
        assert set(restricted) <= {(0, 0, 0, 0)}
        values.append(restricted.get((0, 0, 0, 0), 0))
    return tuple(values)


def check_quartet_unimodularity():
    X = [(1, ((0, 1), (2, 3)))]
    Z = [(1, ((0, 3), (1, 2)))]
    basis = [multiply_terms(X, X), multiply_terms(X, Z), multiply_terms(Z, Z)]
    columns = [quartet_balanced_values(terms) for terms in basis]
    matrix = [[columns[col][row] for col in range(3)] for row in range(3)]
    assert matrix == [[0, 0, 1], [1, -1, 1], [1, 0, 0]]
    determinant = (matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
                   - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
                   + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]))
    assert determinant == 1


def complement_k1_terms(B, variant):
    if variant == 1:
        groups = (B[:3], B[3:6])
    else:
        groups = ((B[0], B[1], B[3]), (B[2], B[4], B[5]))
    terms = multiply_terms(triangle(*groups[0]), triangle(*groups[1]))
    for a, b in zip(B[6::2], B[7::2]):
        terms = multiply_terms(terms, [(1, ((a, b), (a, b)))])
    return terms


def k2_numeric_zero_terms(m, directions):
    assert m == 12
    groups1 = ((0, 1, 2), (3, 4, 5), (6, 7, 8), (9, 10, 11))
    groups2 = ((0, 1, 3), (2, 4, 6), (5, 7, 9), (8, 10, 11))

    def product_of_triangles(groups):
        terms = [(1, ())]
        for group in groups:
            terms = multiply_terms(terms, triangle(*group))
        return terms

    H1 = product_of_triangles(groups1)
    H2 = product_of_triangles(groups2)
    value1, value2 = terms_value(H1, directions), terms_value(H2, directions)
    assert value1 and value2
    relation = add_terms(scale_terms(H1, value2), scale_terms(H2, -value1))
    assert terms_value(relation, directions) == 0
    return relation


def check_coefficient_arrays(m):
    q = m // 2
    B = tuple(range(4, m))
    directions = tuple((0, 1, 2, 3) + tuple(10 + j for j in range(m - 4)))
    Q1, Q2 = quartet_terms()
    assert terms_value(Q1, directions) == terms_value(Q2, directions) == 0
    P1, P2 = complement_k1_terms(B, 1), complement_k1_terms(B, 2)
    c1, c2 = quartet_balanced_values(Q1), quartet_balanced_values(Q2)
    assert c1 == (0, 4, 3) and c2 == (-1, -4, 0)
    W = add_terms(multiply_terms(Q1, P1), scale_terms(multiply_terms(Q2, P2), 2))
    if m == 12:
        W = add_terms(W, k2_numeric_zero_terms(m, directions))
    assert terms_value(W, directions) == 0

    nonzero_arrays = 0
    cut_count = 0
    for S in combinations(range(m), q - 1):
        inside = set(S)
        a = len(inside.intersection(range(4)))
        restricted = restrict_terms(W, inside, m)
        if a != 2:
            assert not restricted, (m, S, a)
        cut_count += 1

    for J in combinations(B, q - 3):
        J = set(J)
        f1 = restrict_terms(P1, J, m)
        f2 = restrict_terms(P2, J, m)
        arrays = [restrict_terms(W, J | {0, i}, m) for i in (1, 2, 3)]
        for i, restricted in zip((1, 2, 3), arrays):
            assert restricted == restrict_terms(W, J | ({0, 1, 2, 3} - {0, i}), m)
        all_monomials = set(f1) | set(f2) | set().union(*(set(a) for a in arrays))
        for monomial in all_monomials:
            p1, p2 = f1.get(monomial, 0), f2.get(monomial, 0)
            array = tuple(a.get(monomial, 0) for a in arrays)
            assert array == tuple(c1[i] * p1 + 2 * c2[i] * p2 for i in range(3))
            if array == (0, 0, 0):
                continue
            nonzero_arrays += 1
            c12, c13, c14 = array
            coeff_X2 = c14
            coeff_XZ = c14 + c12 - c13
            coeff_Z2 = c12
            X_value, Z_value = 1, 3
            assert (coeff_X2 * X_value ** 2
                    + coeff_XZ * X_value * Z_value
                    + coeff_Z2 * Z_value ** 2) == 0

    assert nonzero_arrays > 0
    return cut_count, nonzero_arrays


def gmultiply(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def gconj(a):
    return a[0], -a[1]


def gnorm(a):
    return a[0] * a[0] + a[1] * a[1]


def gpower(a, exponent):
    out = (1, 0)
    for _ in range(exponent):
        out = gmultiply(out, a)
    return out


def nearest_integer(numerator, denominator):
    if numerator < 0:
        return -nearest_integer(-numerator, denominator)
    return (2 * numerator + denominator) // (2 * denominator)


def ggcd(a, b):
    while b != (0, 0):
        denominator = gnorm(b)
        real = nearest_integer(a[0] * b[0] + a[1] * b[1], denominator)
        imag = nearest_integer(a[1] * b[0] - a[0] * b[1], denominator)
        quotient = (real, imag)
        remainder = gmultiply(quotient, b)
        a, b = b, (a[0] - remainder[0], a[1] - remainder[1])
    return a


def gaussian_divides(divisor, value):
    numerator = gmultiply(value, gconj(divisor))
    denominator = gnorm(divisor)
    return all(entry % denominator == 0 for entry in numerator)


def is_prime(p):
    if p < 2:
        return False
    for d in range(2, isqrt(p) + 1):
        if p % d == 0:
            return False
    return True


def split_generators(count):
    out = []
    p = 5
    while len(out) < count:
        if p % 4 == 1 and is_prime(p):
            for a in range(1, isqrt(p) + 1):
                b2 = p - a * a
                b = isqrt(b2)
                if b > 0 and b * b == b2:
                    out.append((a, b))
                    break
            else:
                raise AssertionError(p)
        p += 1
    return out


def bracket(a, b):
    return a[0] * b[1] - a[1] * b[0]


def check_gaussian_aggregation(m):
    A = set(range(4))
    subsets = [frozenset(T) for size in range(1, m + 1)
               for T in combinations(range(m), size)]
    generators = split_generators(len(subsets))
    special = frozenset((1, 2))
    blocks = {}
    for j, (T, pi) in enumerate(zip(subsets, generators)):
        exponent = 3 if T == special else 1 + j % 3
        blocks[T] = gpower(pi, exponent)
    special_pi = generators[subsets.index(special)]
    K = [(1, 0) for _ in range(m)]
    K[0] = special_pi
    K[3] = gpower(special_pi, 2)

    rows = []
    for i in range(m):
        value = K[i]
        for T in subsets:
            if i in T:
                value = gmultiply(value, blocks[T])
        rows.append(value)
        assert gnorm(ggcd(value, gconj(value))) == 1

    X = bracket(rows[0], rows[1]) * bracket(rows[2], rows[3])
    Y = bracket(rows[0], rows[2]) * bracket(rows[1], rows[3])
    Z = bracket(rows[0], rows[3]) * bracket(rows[1], rows[2])
    assert X - Y + Z == 0 and X and Y
    content = gcd(abs(X), abs(Y))
    u, v = X // content, Y // content
    c12, c13, c14 = -u * u, v * v, v * v - u * u
    assert gcd(abs(c12), abs(c13), abs(c14)) == 1
    values = {(0, 1): c12, (0, 2): c13, (0, 3): c14}
    checked = 0
    for I, c_I in values.items():
        paired_moduli = []
        for pair in (set(I), A - set(I)):
            H = (1, 0)
            for T in subsets:
                if T.intersection(A) == pair:
                    H = gmultiply(H, blocks[T])
            assert sum(T.intersection(A) == pair for T in subsets) == 2 ** (m - 4)
            correction = (1, 0)
            for j in A - pair:
                correction = gmultiply(correction, gpower(K[j], 2))
            assert gaussian_divides(H, gmultiply(correction, (c_I, 0)))
            common = ggcd(H, correction)
            modulus = gnorm(H) // gnorm(common)
            assert c_I % modulus == 0
            paired_moduli.append(modulus)
            checked += 1
        assert gcd(*paired_moduli) == 1
        assert c_I % (paired_moduli[0] * paired_moduli[1]) == 0
    return checked


def main():
    check_quartet_unimodularity()
    cuts10, arrays10 = check_coefficient_arrays(10)
    cuts12, arrays12 = check_coefficient_arrays(12)
    gaussian4 = check_gaussian_aggregation(4)
    gaussian6 = check_gaussian_aggregation(6)
    print("PASS: unimodular four-row balanced map and integral arrays;")
    print(f"      m=10: {cuts10} first-smaller cuts, {arrays10} nonzero arrays;")
    print(f"      m=12: {cuts12} first-smaller cuts, {arrays12} nonzero arrays, numeric K2 term;")
    print(f"      Gaussian aggregate divisibility: m=4 ({gaussian4}) and m=6 ({gaussian6}) pair checks.")


if __name__ == "__main__":
    main()
