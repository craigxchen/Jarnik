"""Exact isolated-pole cancellation on actual common-unit Gaussian circles.

The Hensel construction is proved for every e in the accompanying note.
These finite checks do not assert endpoint-scale examples.
"""

from fractions import Fraction
from itertools import combinations
from math import atan, factorial, gcd, prod


def permanent_coefficients(a, modulus, marked=(), derivative_index=None):
    """Permanent in a formal pole variable, with its first derivative."""
    n = len(a)
    values, derivatives = [], []
    for i in range(n):
        value_row, derivative_row = [], []
        for j in range(n):
            if (i, j) in marked:
                value_row.append((1, 1))
                derivative_row.append(0)
            else:
                inv = pow((1 + a[i] * a[j]) % modulus, -1, modulus)
                derivative = ((a[j] if i == derivative_index else 0)
                              + (a[i] if j == derivative_index else 0))
                value_row.append((0, inv))
                derivative_row.append(-derivative * inv * inv % modulus)
        values.append(value_row)
        derivatives.append(derivative_row)
    size = 1 << n
    dp = [[0] * 3 for _ in range(size)]
    dd = [[0] * 3 for _ in range(size)]
    dp[0][0] = 1
    for mask in range(size):
        row = bin(mask).count('1')
        if row == n:
            continue
        for col in range(n):
            if mask >> col & 1:
                continue
            nxt = mask | (1 << col)
            shift, value = values[row][col]
            deriv = derivatives[row][col]
            for degree in range(3 - shift):
                target = degree + shift
                dp[nxt][target] = (dp[nxt][target] + dp[mask][degree] * value) % modulus
                dd[nxt][target] = (dd[nxt][target] + dd[mask][degree] * value
                                  + dp[mask][degree] * deriv) % modulus
    return dp[-1], dd[-1]


def F_value(a, modulus):
    delta = 1 + a[0] * a[1]
    coeff, deriv = permanent_coefficients(a, modulus, ((0, 1), (1, 0)), 7)
    value = (delta * delta * coeff[0] + delta * coeff[1] + coeff[2]) % modulus
    slope = (delta * delta * deriv[0] + delta * deriv[1] + deriv[2]) % modulus
    return value, slope


def crt(r, modulus, s, other):
    assert gcd(modulus, other) == 1
    return r + modulus * (((s - r) * pow(modulus, -1, other)) % other)


def vp(integer, prime):
    assert integer
    answer = 0
    while integer % prime == 0:
        integer //= prime
        answer += 1
    return answer


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def mul(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def conj(a):
    return a[0], -a[1]


def norm(a):
    return a[0] ** 2 + a[1] ** 2


def divide(a, b):
    numerator = mul(a, conj(b))
    denominator = norm(b)
    assert all(x % denominator == 0 for x in numerator)
    return tuple(x // denominator for x in numerator)


def ggcd(a, b):
    while b != (0, 0):
        numerator = mul(a, conj(b))
        denominator = norm(b)
        quotient = tuple((2 * x + denominator) // (2 * denominator) for x in numerator)
        remainder = add(a, tuple(-x for x in mul(quotient, b)))
        assert norm(remainder) < norm(b)
        a, b = b, remainder
    return a


def is_primary(a):
    return a[1] % 2 == 0 and (a[0] + a[1]) % 4 == 1


def primary_associate(a):
    candidates = [mul(a, unit) for unit in ((1, 0), (-1, 0), (0, 1), (0, -1))]
    primary = [value for value in candidates if is_primary(value)]
    assert len(primary) == 1
    return primary[0]


base = [12, 2, 6, 11, 15, 5]
coeff, derivative = permanent_coefficients(base, 17, derivative_index=5)
assert coeff[0] == 0 and derivative[0] == 2
for e in range(1, 5):
    pole_modulus = 17 ** (2 * e + 1)
    pole_residue = (17 ** (2 * e) - 1) * pow(52, -1, pole_modulus) % pole_modulus
    a2 = crt(pole_residue, pole_modulus, 0, 4)
    a = [52, a2, 12, 36, 40, 28, 32, 5]
    assert vp(1 + 52 * a2, 17) == 2 * e
    for k in range(1, 4 * e):
        modulus = 17 ** (k + 1)
        value, derivative = F_value(a, modulus)
        assert value % (17 ** k) == 0 and derivative % 17 == 2
        digit = -(value // (17 ** k)) * pow(derivative % 17, -1, 17) % 17
        a[7] += digit * (17 ** k)
    a[7] = crt(a[7], 17 ** (4 * e), 0, 4)
    assert F_value(a, 17 ** (4 * e))[0] == 0
    assert all(value > 0 and value % 4 == 0 for value in a)
    assert [value % 17 for value in a] == [1, 16, 12, 2, 6, 11, 15, 5]
    assert all((value * value + 1) % 17 for value in a)
    assert all((a[i] - a[j]) % 17 for i, j in combinations(range(8), 2))
    assert [(i, j) for i, j in combinations(range(8), 2)
            if (1 + a[i] * a[j]) % 17 == 0] == [(0, 1)]

    h = [(value, 1) for value in a]
    assert all(norm(ggcd(value, conj(value))) == 1 for value in h)
    for value in h:
        primary = mul((0, -1), value)
        divide(add(primary, (-1, 0)), (-2, 2))  # (1+i)^3=-2+2i.
        assert is_primary(primary)
    L = (1, 0)
    for value in h:
        L = mul(divide(L, ggcd(L, value)), value)
    L = primary_associate(L)
    points = [divide(mul(conj(L), value), conj(value)) for value in h]
    assert len(set(points)) == 8
    assert all(norm(point) == norm(L) for point in points)
    common = (0, 0)
    for point in points:
        common = ggcd(common, point)
    common = primary_associate(common)
    points = [divide(point, common) for point in points]
    assert all(is_primary(tuple(-x for x in point)) for point in points)
    common = (0, 0)
    for point in points:
        common = ggcd(common, point)
    assert norm(common) == 1
    assert len({norm(point) for point in points}) == 1
    assert norm(points[0]) % 2 and norm(points[0]) % 17
    for i, j in combinations(range(8), 2):
        overlap = norm(ggcd(h[i], h[j]))
        relative = divide(mul(h[i], conj(h[j])), (overlap, 0))
        xx, tt = abs(relative[0]), abs(relative[1])
        assert gcd(xx, tt) == 1 and tt % 17
        assert vp(xx, 17) == (2 * e if (i, j) == (0, 1) else 0)
    assert max(2 * atan(1 / value) for value in a) < 1

    if e == 1:
        # One full reduced rational permanent, independently of modular evaluation.
        dp = [Fraction(0)] * 256
        dp[0] = Fraction(1)
        for mask in range(256):
            i = bin(mask).count('1')
            if i == 8:
                continue
            for j in range(8):
                if not(mask >> j & 1):
                    dp[mask | (1 << j)] += dp[mask] / (1 + a[i] * a[j])
        permanent = dp[-1] * prod(value * value + 1 for value in a)
        assert permanent.denominator % 17 and permanent > factorial(8)
        first_taylor = factorial(8) + factorial(7) * sum(
            Fraction((a[i] - a[j]) ** 2, (a[i] ** 2 + 1) * (a[j] ** 2 + 1))
            for i, j in combinations(range(8), 2))
        remainder = permanent - first_taylor
        assert remainder > 0 and remainder.denominator % 17

print("Complementary six-row permanent: value0 and derivative2 modulo17 verified.")
print("Hensel construction checked for e=1,2,3,4, with exact pole orders2e.")
print("Actual primitive odd-norm Gaussian circles and every old residue verified.")
print("One fixed primary-unit convention and full local denominator cancellation pass.")
print("One full rational permanent and Taylor remainder independently reduced exactly.")
print("These fixed-span families are not endpoint-scale counterexamples.")
