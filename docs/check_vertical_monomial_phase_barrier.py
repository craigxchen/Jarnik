"""Exact integer checks for vertical_monomial_phase_barrier.md."""

from itertools import combinations, product
from math import gcd

from check_endpoint_vertical_line_gcd import exact_div, gcd_gaussian, norm


def mul(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def power(a, n):
    result = (1, 0)
    for _ in range(n):
        result = mul(result, a)
    return result


def check_signed_incidence():
    subsets = [s for size in (2, 3) for s in combinations(range(4), size)]
    count = 0
    minimum = None
    for r in product(range(-4, 5), repeat=4):
        if not any(r):
            continue
        total = sum(abs(sum(r[j] for j in subset)) for subset in subsets)
        assert total >= 6
        minimum = total if minimum is None else min(minimum, total)
        count += 1
    assert minimum == 6
    return count


def check_quotient(d, ts, r):
    u = v = (1, 0)
    for t, exponent in zip(ts, r):
        if exponent > 0:
            u = mul(u, power((-d, t), exponent))
        elif exponent < 0:
            v = mul(v, power((-d, t), -exponent))
    common = gcd_gaussian(u, v)
    n, denominator = exact_div(u, common), exact_div(v, common)
    z = mul(n, (denominator[0], -denominator[1]))
    assert norm(z) * norm(common)**2 == norm(u) * norm(v)
    assert norm(gcd_gaussian(n, denominator)) == 1
    # Multiplication by i^(-sum r), without any floating-point angles.
    rotated = mul(z, power((0, 1), (-sum(r)) % 4))
    content = gcd(abs(rotated[0]), abs(rotated[1]))
    primitive = (rotated[0] // content, rotated[1] // content)
    weight = sum(map(abs, r))
    for value in (rotated, primitive):
        assert value[1]**2 * min(ts)**2 <= norm(value) * d**2 * weight**2
        if value[1]:
            assert norm(value) * d**2 * weight**2 >= min(ts)**2
    return rotated[1] == 0


def check_actual_quotients():
    count = zero_count = 0
    for d in (1, 2, 5):
        for ts in ((1, 2, 3, 4), (7, 11, 19, 31), (100, 103, 109, 127)):
            for r in product(range(-2, 3), repeat=4):
                if not any(r):
                    continue
                zero_count += check_quotient(d, ts, r)
                count += 1
    # This deliberately preserves the zero-projection exception.
    # (-1+i)(-1+2i)(-1+3i)=10, whose square has s=6.
    assert check_quotient(1, (1, 2, 3, 4), (2, 2, 2, 0))
    assert zero_count > 0
    return count, zero_count


if __name__ == "__main__":
    print(f"PASS: {check_signed_incidence()} signed incidence vectors; minimum sum 6")
    count, zero_count = check_actual_quotients()
    print(f"PASS: {count} exact Gaussian quotient checks, including {zero_count} zero projections")
