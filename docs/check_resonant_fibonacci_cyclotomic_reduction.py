"""Exact small fixtures for resonant_fibonacci_cyclotomic_reduction.md."""

from fractions import Fraction
from math import gcd


def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def lucas(n):
    return 2 if n == 0 else fib(n - 1) + fib(n + 1)


def gaussian(d):
    assert d > 0 and d % 2 == 1
    return fib((d + 1) // 2), fib((d - 1) // 2)


def norm(z):
    return z[0] * z[0] + z[1] * z[1]


def conjugate(z):
    return z[0], -z[1]


def multiply(z, w):
    a, b = z
    c, d = w
    return a * c - b * d, a * d + b * c


def gaussian_power(z, exponent):
    answer = (1, 0)
    for _ in range(exponent):
        answer = multiply(answer, z)
    return answer


def quotient_remainder(z, w):
    a, b = z
    c, d = w
    divisor_norm = norm(w)
    q = (
        round(Fraction(a * c + b * d, divisor_norm)),
        round(Fraction(b * c - a * d, divisor_norm)),
    )
    qw = multiply(q, w)
    return q, (a - qw[0], b - qw[1])


def exact_quotient(z, w):
    q, remainder = quotient_remainder(z, w)
    assert remainder == (0, 0)
    return q


def gaussian_gcd(z, w):
    while w != (0, 0):
        z, w = w, quotient_remainder(z, w)[1]
    return z


def polynomial_multiply(a, b):
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] += x * y
    return result


def polynomial_exact_quotient(numerator, denominator):
    result = [0] * (len(numerator) - len(denominator) + 1)
    remainder = numerator[:]
    for position in range(len(result) - 1, -1, -1):
        coefficient = remainder[position + len(denominator) - 1]
        result[position] = coefficient
        for j, value in enumerate(denominator):
            remainder[position + j] -= coefficient * value
    assert all(value == 0 for value in remainder)
    return result


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def cyclotomic(n, cache):
    if n in cache:
        return cache[n]
    polynomial = [-1] + [0] * (n - 1) + [1]
    for d in divisors(n)[:-1]:
        polynomial = polynomial_exact_quotient(polynomial, cyclotomic(d, cache))
    cache[n] = polynomial
    return polynomial


def power_i(exponent):
    return ((1, 0), (0, 1), (-1, 0), (0, -1))[exponent % 4]


def homogeneous_layer(e, d0, cache):
    """Evaluate Phi_e(alpha^d0,beta^d0) using Gaussian Lucas sums."""
    assert e > 1 and e % 2 == 1 and d0 % 2 == 1
    coefficients = cyclotomic(e, cache)
    degree = len(coefficients) - 1
    assert degree % 2 == 0 and coefficients == coefficients[::-1]
    base_product = power_i(d0)
    answer = (0, 0)
    for j in range(degree // 2):
        term = power_i(d0 * j)
        scalar = coefficients[j] * lucas((degree // 2 - j) * d0)
        answer = answer[0] + scalar * term[0], answer[1] + scalar * term[1]
    middle = power_i(d0 * (degree // 2))
    answer = (
        answer[0] + coefficients[degree // 2] * middle[0],
        answer[1] + coefficients[degree // 2] * middle[1],
    )
    assert norm(base_product) == 1
    return answer


def mobius(n):
    count = 0
    p = 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0:
                return 0
            count += 1
        p += 1
    if n > 1:
        count += 1
    return -1 if count % 2 else 1


def ramanujan(e, k):
    return sum(d * mobius(e // d) for d in divisors(gcd(e, k)))


def check_strong_divisibility():
    count = 0
    for d in range(1, 70, 2):
        assert norm(gaussian(d)) == fib(d)
        assert norm(gaussian_gcd(gaussian(d), conjugate(gaussian(d)))) <= 4
        for e in range(1, 70, 2):
            actual = norm(gaussian_gcd(gaussian(d), gaussian(e)))
            assert actual == fib(gcd(d, e))
            count += 1
    for n in range(1, 7):
        d = 12 * n + 1
        q = exact_quotient(gaussian(13 * d), gaussian(d))
        assert norm(q) == fib(13 * d) // fib(d)
    print(f"PASS: {count} oriented Gaussian gcd pairs and six 13:1 quotients")


def check_layers():
    cache = {}
    tested = 0
    for d0 in (13, 25, 37):
        layers = {1: gaussian(d0)}
        for e in (3, 5, 7, 9, 15):
            layer = homogeneous_layer(e, d0, cache)
            lower = (1, 0)
            for f in divisors(e)[:-1]:
                lower = multiply(lower, layers[f])
            assert multiply(lower, layer) == gaussian(e * d0)
            layers[e] = layer
            tested += 1
        for e in layers:
            for f in layers:
                if e < f:
                    assert norm(gaussian_gcd(layers[e], layers[f])) <= 5
                assert norm(gaussian_gcd(layers[e], conjugate(layers[f]))) <= 4
    print(f"PASS: {tested} integral cyclotomic layers and bounded-overlap fixtures")


def check_effective_radius():
    cache = {}
    multiplicities = (1, 1, 1)
    multiples = (1, 3, 9)
    rows = ((0, 0, 0), (1, 0, 1), (0, 1, 1))
    layer_indices = (1, 3, 9)
    for d0 in (13, 25, 37):
        layers = {1: gaussian(d0)}
        for e in layer_indices[1:]:
            layers[e] = homogeneous_layer(e, d0, cache)
        exponents = {
            e: [sum(row[s] for s, m in enumerate(multiples) if m % e == 0)
                for row in rows]
            for e in layer_indices
        }
        totals = {
            e: sum(multiplicities[s] for s, m in enumerate(multiples) if m % e == 0)
            for e in layer_indices
        }
        common = (1, 0)
        predicted_norm = 1
        for e in layer_indices:
            values = exponents[e]
            width = max(values) - min(values)
            predicted_norm *= norm(layers[e]) ** width
            for _ in range(min(values)):
                common = multiply(common, layers[e])
            for _ in range(totals[e] - max(values)):
                common = multiply(common, conjugate(layers[e]))
        actual_rows = []
        for row in rows:
            value = (1, 0)
            for s, m in enumerate(multiples):
                block = gaussian(m * d0)
                for _ in range(row[s]):
                    value = multiply(value, block)
                for _ in range(multiplicities[s] - row[s]):
                    value = multiply(value, conjugate(block))
            actual_rows.append(value)
        quotients = [exact_quotient(value, common) for value in actual_rows]
        assert len({norm(value) for value in quotients}) == 1
        assert norm(quotients[0]) == predicted_norm
        whole_gcd = actual_rows[0]
        for value in actual_rows[1:]:
            whole_gcd = gaussian_gcd(whole_gcd, value)
        excess = exact_quotient(whole_gcd, common)
        assert norm(excess) <= 4
    print("PASS: explicit common layer factor and primitive-radius fixtures")


def check_polynomial_fibre():
    cache = {}
    for r in (1, 3, 5, 7):
        plus = [0] * (2 * r + 1)
        plus[0] = plus[r] = plus[2 * r] = 1
        minus = plus[:]
        minus[r] = -1
        product = [1]
        for e in divisors(3 * r):
            if r % e != 0 and e != 1:
                product = polynomial_multiply(product, cyclotomic(e, cache))
        assert product == plus
        assert [x * ((-1) ** j) for j, x in enumerate(product)] == minus
        assert polynomial_multiply(plus, [x * ((-1) ** j) for j, x in enumerate(plus)]) == polynomial_multiply(minus, [x * ((-1) ** j) for j, x in enumerate(minus)])
        assert next(j for j, x in enumerate(a - b for a, b in zip(plus, minus)) if x) == r
        for d0 in (1, 13, 25):
            k = r * d0
            quotient = exact_quotient(gaussian(3 * k), gaussian(k))
            assert quotient == (
                lucas(k),
                power_i(k)[1],
            )
            pair_gcd = gaussian_gcd(quotient, conjugate(quotient))
            assert norm(pair_gcd) == (1 if k % 3 == 0 else 2)

    b = {3: 1, 15: 1}
    for k in (1, 3):
        assert sum(value * ramanujan(e, k) for e, value in b.items()) == 0
    assert sum(value * ramanujan(e, 5) for e, value in b.items()) != 0
    for d in (1, 3):
        assert sum(value * mobius(e // d) for e, value in b.items() if e % d == 0) == 0
    print("PASS: reciprocal fibre, contact orders, and Ramanujan inversion")


def check_converse_incidence():
    """Invert arbitrary layer orientations and lift them to raw Gaussian rows."""
    cache = {}
    layer_indices = (1, 3, 5, 15)
    prescribed = (
        {1: 0, 3: 0, 5: 0, 15: 0},
        {1: 2, 3: 1, 5: 0, 15: 1},
        {1: 0, 3: 0, 5: 1, 15: 1},
    )
    raw = [
        {
            m: sum(mobius(e // m) * row[e]
                   for e in layer_indices if e % m == 0)
            for m in layer_indices
        }
        for row in prescribed
    ]
    shifts = {m: -min(row[m] for row in raw) for m in layer_indices}
    multiplicities = {
        m: max(row[m] for row in raw) - min(row[m] for row in raw)
        for m in layer_indices
    }
    orientations = [
        {m: row[m] + shifts[m] for m in layer_indices}
        for row in raw
    ]
    for i, row in enumerate(orientations):
        assert all(0 <= row[m] <= multiplicities[m] for m in layer_indices)
        for e in layer_indices:
            recovered = sum(row[m] for m in layer_indices if m % e == 0)
            common_offset = sum(shifts[m] for m in layer_indices if m % e == 0)
            assert recovered == prescribed[i][e] + common_offset

    totals = {e: sum(multiplicities[m] for m in layer_indices if m % e == 0)
              for e in layer_indices}
    layer_counts = {
        e: [sum(row[m] for m in layer_indices if m % e == 0)
            for row in orientations]
        for e in layer_indices
    }
    base_counts = layer_counts[1]
    assert all((value - base_counts[0]) % 2 == 0 for value in base_counts)
    half_differences = [(value - base_counts[0]) // 2 for value in base_counts]
    K = max(abs(k) for k in half_differences)
    F, Fbar = (1, 2), (1, -2)
    prefactors = [
        multiply(gaussian_power(F, K - k), gaussian_power(Fbar, K + k))
        for k in half_differences
    ]
    assert len({norm(value) for value in prefactors}) == 1
    for i, k in enumerate(half_differences):
        assert multiply(prefactors[i], gaussian_power(F, k)) == (
            multiply(prefactors[0], gaussian_power(Fbar, k))
        )

    for d0 in (13, 25):
        layers = {1: gaussian(d0)}
        for e in layer_indices[1:]:
            layers[e] = homogeneous_layer(e, d0, cache)
        common = (1, 0)
        for e in layer_indices:
            counts = layer_counts[e]
            common = multiply(
                common,
                gaussian_power(layers[e], min(counts)),
            )
            common = multiply(
                common,
                gaussian_power(conjugate(layers[e]), totals[e] - max(counts)),
            )
        normalized_norms = []
        for i, row in enumerate(orientations):
            actual = prefactors[i]
            for m in layer_indices:
                actual = multiply(
                    actual,
                    gaussian_power(gaussian(m * d0), row[m]),
                )
                actual = multiply(
                    actual,
                    gaussian_power(conjugate(gaussian(m * d0)),
                                   multiplicities[m] - row[m]),
                )
            normalized = exact_quotient(actual, common)
            expected = prefactors[i]
            for e in layer_indices:
                count = layer_counts[e][i] - min(layer_counts[e])
                width = max(layer_counts[e]) - min(layer_counts[e])
                expected = multiply(expected, gaussian_power(layers[e], count))
                expected = multiply(
                    expected,
                    gaussian_power(conjugate(layers[e]), width - count),
                )
            assert normalized == expected
            normalized_norms.append(norm(normalized))
        assert len(set(normalized_norms)) == 1
    print("PASS: arbitrary layer incidence inversion and Gaussian lift")


if __name__ == "__main__":
    check_strong_divisibility()
    check_layers()
    check_effective_radius()
    check_polynomial_fibre()
    check_converse_incidence()
