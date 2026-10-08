"""Exact finite arithmetic checks for actual-cut nonvanishing.

The m=6 terminal and m=8,k=1 nonterminal fixtures use a distinct odd
split prime for every nonempty proper full-core block. They are actual
Gaussian rows, not short-arc examples and not near-uniform profiles.
The checker verifies the projected-block divisibility mechanism; the
all-depth theorem is proved symbolically in its accompanying note.
"""

from itertools import combinations
from math import gcd, isqrt


def gmul(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def gconj(a):
    return a[0], -a[1]


def gnorm(a):
    return a[0] * a[0] + a[1] * a[1]


def gpow(a, power):
    out = (1, 0)
    for _ in range(power):
        out = gmul(out, a)
    return out


def nearest(numerator, denominator):
    if numerator < 0:
        return -nearest(-numerator, denominator)
    return (2 * numerator + denominator) // (2 * denominator)


def ggcd(a, b):
    while b != (0, 0):
        denominator = gnorm(b)
        real = nearest(a[0] * b[0] + a[1] * b[1], denominator)
        imag = nearest(a[1] * b[0] - a[0] * b[1], denominator)
        product = gmul((real, imag), b)
        a, b = b, (a[0] - product[0], a[1] - product[1])
    return a


def prime(p):
    if p < 2:
        return False
    return all(p % d for d in range(2, isqrt(p) + 1))


def split_generators(count):
    answer = []
    p = 5
    while len(answer) < count:
        if p % 4 == 1 and prime(p):
            for a in range(1, isqrt(p) + 1):
                b2 = p - a * a
                b = isqrt(b2)
                if b and b * b == b2:
                    answer.append((a, b))
                    break
            else:
                raise AssertionError(p)
        p += 1
    return answer


def core_fixture(m, inside, chosen_y_set):
    outside = tuple(i for i in range(m) if i not in inside)
    subsets = [frozenset(T) for size in range(1, m)
               for T in combinations(range(m), size)]
    primes = split_generators(len(subsets))
    special = frozenset(chosen_y_set)
    assert special in subsets
    pi_special = primes[subsets.index(special)]
    blocks = {T: gpow(pi, 3 if T == special else 1)
              for T, pi in zip(subsets, primes)}

    correction = [(1, 0) for _ in range(m)]
    for i in tuple(j for j in outside if j not in chosen_y_set)[:2]:
        correction[i] = pi_special
    rows = []
    for i in range(m):
        row = correction[i]
        for T, block in blocks.items():
            if i in T:
                row = gmul(row, block)
        rows.append(row)
        assert gnorm(ggcd(row, gconj(row))) == 1
    return outside, blocks, correction, rows, special, pi_special


def determinant(rows, a, b):
    # P_a Y_b-P_b Y_a is real: its imaginary terms cancel.
    value = rows[a][0] * rows[b][1] - rows[b][0] * rows[a][1]
    assert value
    return value


def covariant_value(covariant, rows):
    unpaired, edges = covariant
    value = (1, 0) if unpaired is None else rows[unpaired]
    for a, b in edges:
        value = gmul(value, (determinant(rows, a, b), 0))
    return value


def covariant_coefficients(covariant):
    # Coefficient of Y_J P_(outside\J), represented by J.
    unpaired, edges = covariant
    coefficients = {frozenset(): 1}
    for a, b in edges:
        next_coefficients = {}
        for J, c in coefficients.items():
            next_coefficients[J | {b}] = next_coefficients.get(J | {b}, 0) + c
            next_coefficients[J | {a}] = next_coefficients.get(J | {a}, 0) - c
        coefficients = next_coefficients
    assert all(len(J) == len(edges) for J in coefficients)
    if unpaired is not None:
        assert all(unpaired not in J for J in coefficients)
    return {J: c for J, c in coefficients.items() if c}


def linear_combination(values, coefficients):
    return (sum(c * value[0] for c, value in zip(coefficients, values)),
            sum(c * value[1] for c, value in zip(coefficients, values)))


def coefficient_dictionary(covariants, coefficients):
    answer = {}
    for covariant, scalar in zip(covariants, coefficients):
        for J, coefficient in covariant_coefficients(covariant).items():
            answer[J] = answer.get(J, 0) + scalar * coefficient
    return {J: coefficient for J, coefficient in answer.items() if coefficient}


def isolated_monomial_value(covariants, coefficients, J):
    # Set (P_i,Y_i)=(0,1) for i in J and (1,0) otherwise.
    value = 0
    for scalar, (unpaired, edges) in zip(coefficients, covariants):
        term = 1
        if unpaired is not None and unpaired in J:
            term = 0
        for a, b in edges:
            p_a, y_a = (0, 1) if a in J else (1, 0)
            p_b, y_b = (0, 1) if b in J else (1, 0)
            term *= p_a * y_b - p_b * y_a
        value += scalar * term
    return value


def projected_block(blocks, outside, J, expected_count):
    selected = [T for T in blocks if T.intersection(outside) == J]
    assert len(selected) == expected_count
    product = (1, 0)
    for T in selected:
        product = gmul(product, blocks[T])
    return product


def verify_residual(m, inside, covariants, coefficients):
    assert m in (6, 8)
    outside = tuple(i for i in range(m) if i not in inside)
    selected_y = frozenset(outside[:2])
    outside, blocks, correction, rows, special, pi_special = core_fixture(
        m, inside, selected_y)
    assert len(blocks) == 2 ** m - 2
    assert all(determinant(rows, i, j)
               for i, j in combinations(range(m), 2))
    values = [covariant_value(covariant, rows) for covariant in covariants]
    if m == 8:
        assert coefficients is None and len(values) == 3
        real = [value[0] for value in values]
        imag = [value[1] for value in values]
        coefficients = [real[1] * imag[2] - imag[1] * real[2],
                        real[2] * imag[0] - imag[2] * real[0],
                        real[0] * imag[1] - imag[0] * real[1]]
    else:
        assert coefficients is None and m == 6 and len(values) == 2
        assert values[0][1] == values[1][1] == 0
        coefficients = [values[1][0], -values[0][0]]
    content = gcd(*[abs(c) for c in coefficients])
    assert content
    coefficients = [c // content for c in coefficients]
    assert linear_combination(values, coefficients) == (0, 0)

    polynomial = coefficient_dictionary(covariants, coefficients)
    assert polynomial
    assert all(len(J) == 2 for J in polynomial)
    for J_tuple in combinations(outside, 2):
        J = frozenset(J_tuple)
        c_J = polynomial.get(J, 0)
        assert isolated_monomial_value(covariants, coefficients, J) == c_J
        H = projected_block(blocks, outside, J, 2 ** len(inside))
        if not c_J:
            continue
        kappa = (1, 0)
        for i in outside:
            if i not in J:
                kappa = gmul(kappa, correction[i])
        common = ggcd(H, kappa)
        modulus = gnorm(H) // gnorm(common)
        assert c_J % modulus == 0
        if m == 6:
            complement = frozenset(outside) - J
            assert polynomial[complement] == c_J
            H_complement = projected_block(blocks, outside, complement,
                                           2 ** len(inside))
            kappa_complement = (1, 0)
            for i in outside:
                if i not in complement:
                    kappa_complement = gmul(kappa_complement, correction[i])
            common_complement = ggcd(H_complement, kappa_complement)
            modulus_complement = (gnorm(H_complement)
                                  // gnorm(common_complement))
            assert gcd(modulus, modulus_complement) == 1
            assert c_J % (modulus * modulus_complement) == 0

    special_H = projected_block(blocks, outside, selected_y, 2 ** len(inside))
    special_kappa = (1, 0)
    for i in outside:
        if i not in selected_y:
            special_kappa = gmul(special_kappa, correction[i])
    special_loss = gnorm(ggcd(special_H, special_kappa))
    assert special_loss == gnorm(pi_special) ** 2
    return len(blocks), 2 ** len(inside), len(polynomial), len(str(max(abs(c) for c in polynomial.values())))


def main():
    terminal = ((None, ((2, 3), (4, 5))),
                (None, ((2, 4), (3, 5))))
    nonterminal = ((7, ((3, 4), (5, 6))),
                   (6, ((3, 5), (4, 7))),
                   (5, ((3, 6), (4, 7))))
    m6 = verify_residual(6, {0, 1}, terminal, None)
    m8 = verify_residual(8, {0, 1, 2}, nonterminal, None)
    print("PASS: projected Gaussian core isolation and exact numeric-zero residuals")
    print(f"      m=6 terminal: {m6[0]} blocks, {m6[1]} per projection,"
          f" {m6[2]} nonzero coefficients ({m6[3]} digits max)")
    print(f"      m=8,k=1: {m8[0]} blocks, {m8[1]} per projection,"
          f" {m8[2]} nonzero coefficients ({m8[3]} digits max)")
    print("      shared-core corrections and terminal complementary norm products verified")


if __name__ == "__main__":
    main()
