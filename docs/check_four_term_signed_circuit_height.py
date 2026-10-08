"""Exact polynomial, Gaussian-divisor and valuation checks for signed circuits."""

from itertools import combinations, product
from math import gcd

from check_norm_descent_rational_phase_gap import add, conj, gaussian_gcd, mul, neg, norm
from check_signed_block_flip_witness_height import exact_quotient, gproduct


def polynomial_product(left, right):
    out = [(0, 0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i+j] = add(out[i+j], mul(a, b))
    return out


def polynomial_conjugate(poly):
    return [conj(a) for a in poly]


def check_polynomial_identity():
    a = [(0, 5), (1, 1)]
    b = [(1, 4), (1, 1)]
    c = [(2, 3), (1, 1)]
    signed = []
    for bits in ((0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0)):
        poly = [(1, 0)]
        for factor, bit in zip((a, b, c), bits):
            poly = polynomial_product(poly, polynomial_conjugate(factor) if bit else factor)
        signed.append(poly)
    expected = [
        [(-55, -50), (-71, -19), (-24, 6), (-2, 2)],
        [(55, -50), (51, -41), (16, -14), (2, -2)],
        [(25, -70), (29, -59), (12, -18), (2, -2)],
        [(-25, -70), (-1, -69), (8, -22), (2, -2)],
    ]
    assert signed == expected
    for degree in range(4):
        value = (0, 0)
        for coefficient, poly in zip((-1, -6, 8, -3), signed):
            value = add(value, mul((coefficient, 0), poly[degree]))
        assert value == (0, 0)


def relation_value(coefficients, values):
    result = (0, 0)
    for coefficient, value in zip(coefficients, values):
        result = add(result, mul(coefficient, value))
    return result


def determinant(z, w):
    return mul(conj(z), w)[1]


def check_families_and_minors():
    divisor_checks = 0
    for index in range(64):
        t = 1 + 15*index
        a, b, c = (t, t+5), (t+1, t+4), (t+2, t+3)
        blocks = (a, b, c)
        norms = tuple(norm(z) for z in blocks)
        assert norms == (2*t*t+10*t+25, 2*t*t+10*t+17, 2*t*t+10*t+13)
        for z in blocks:
            assert norm(z) > 1 and norm(gaussian_gcd(z, conj(z))) == 1
        assert all(gcd(x, y) == 1 for x, y in combinations(norms, 2))
        values = [gproduct((a, b, c)), gproduct((a, conj(b), conj(c))),
                  gproduct((conj(a), b, conj(c))), gproduct((conj(a), conj(b), c))]
        assert len(set(values)) == 4
        assert all(norm(z) == norms[0]*norms[1]*norms[2] for z in values)
        common = values[0]
        for z in values[1:]:
            common = gaussian_gcd(common, z)
        assert norm(common) == 1
        aligned = [neg(values[0])] + values[1:]
        assert all(determinant(aligned[i], aligned[j]) < 0
                   for i, j in combinations(range(4), 2))

        first = [(-1, 0), (-6, 0), (8, 0), (-3, 0)]
        second = [(0, 0), neg(conj(a)), mul((2, 0), conj(b)), neg(conj(c))]
        real_second = [(determinant(values[1], values[2]), 0),
                       (-determinant(values[0], values[2]), 0),
                       (determinant(values[0], values[1]), 0), (0, 0)]
        assert add(add(a, mul((-2, 0), b)), c) == (0, 0)
        assert relation_value(first, values) == (0, 0)
        divisors = {(0, 1): conj(a), (2, 3): a, (0, 2): conj(b),
                    (1, 3): b, (0, 3): conj(c), (1, 2): c}
        for coefficients, real in ((second, False), (real_second, True)):
            assert relation_value(coefficients, values) == (0, 0)
            nonzero = False
            for (i, j), divisor in divisors.items():
                minor = add(mul(first[i], coefficients[j]), neg(mul(first[j], coefficients[i])))
                exact_quotient(minor, divisor)
                nonzero |= minor != (0, 0)
                if real:
                    assert minor[1] == 0 and minor[0] % norm(divisor) == 0
                divisor_checks += 1
            assert nonzero
            first_height_sq = max(map(norm, first))
            second_height_sq = max(map(norm, coefficients))
            assert 4*first_height_sq*second_height_sq >= min(norms)
            if real:
                assert 2*max(abs(z[0]) for z in first)*max(abs(z[0]) for z in coefficients) >= min(norms)
        assert max(map(norm, second)) == 4*norm(b)
    return divisor_checks


def check_normalization_cost():
    count = 0
    for depth in range(1, 7):
        original = [(depth, 0), (0, depth), (0, 0), (0, 0)]
        for a in range(9):
            for b in range(9):
                for choices in product((0, 1), repeat=4):
                    residues = []
                    for (u, v), choice in zip(original, choices):
                        residues.append((u+a-depth*(choice == 0), v+b-depth*(choice == 1)))
                    if any(min(pair) < 0 for pair in residues):
                        continue
                    total = sum(u+v for u, v in residues)
                    assert a+b >= depth
                    assert total == 4*(a+b)-2*depth and total >= 2*depth
                    count += 1
    return count


if __name__ == '__main__':
    check_polynomial_identity()
    divisors = check_families_and_minors()
    valuations = check_normalization_cost()
    print('PASS: complete cubic polynomial identity and 64 primitive disjoint-support fixtures.')
    print(f'PASS: both independent relations, {divisors} oriented coefficient-minor checks, and height bounds.')
    print(f'PASS: {valuations} admissible prime-depth normalization cases.')
