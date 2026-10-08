"""Exact finite audit of pair-reflection domains on Gaussian circle arcs."""

from functools import cmp_to_key
from itertools import combinations
from math import isqrt

from check_gaussian_reflection_replacement import (
    gconj,
    gdivexact,
    ggcd,
    gmul,
    gnorm,
)


MAX_N = 500


def divisible(divisor, value):
    denominator = gnorm(divisor)
    numerator = gmul(value, gconj(divisor))
    return all(coordinate % denominator == 0 for coordinate in numerator)


def angular_compare(left, right):
    """Exact counterclockwise order starting at the positive real axis."""
    left_half = 0 if left[1] > 0 or (left[1] == 0 and left[0] >= 0) else 1
    right_half = 0 if right[1] > 0 or (right[1] == 0 and right[0] >= 0) else 1
    if left_half != right_half:
        return -1 if left_half < right_half else 1
    cross = left[0] * right[1] - left[1] * right[0]
    if cross:
        return -1 if cross > 0 else 1
    return 0


def circle_points(n):
    bound = isqrt(n)
    points = {
        (x, y)
        for x in range(-bound, bound + 1)
        for y in range(-bound, bound + 1)
        if x * x + y * y == n
    }
    return sorted(points, key=cmp_to_key(angular_compare))


def cyclic_short_arcs(points):
    """All cyclic consecutive blocks whose endpoint separation is < pi/2."""
    count = len(points)
    for start in range(count):
        for length in range(2, count + 1):
            block = tuple(points[(start + offset) % count] for offset in range(length))
            dot = block[0][0] * block[-1][0] + block[0][1] * block[-1][1]
            cross = block[0][0] * block[-1][1] - block[0][1] * block[-1][0]
            if dot <= 0 or cross <= 0:
                break
            yield block


def check_exact_capacity(n, endpoint_dot, D, domain_size, internal_size):
    """Check exact sagitta forms which imply the two Delta bounds."""
    # Keep the ceiling from the number of occupied projection levels.
    # 2(ceil(|S|/2)-1) <= 2 sqrt(N/D) (1-cos Delta).
    domain_excess = 2 * ((domain_size + 1) // 2 - 1)
    assert domain_excess >= 0
    assert domain_excess**2 * n * D <= 4 * (n - endpoint_dot) ** 2

    # 2(ceil(|A|/2)-1) <= 2 sqrt(N/D) (1-cos(Delta/2)). Isolating the
    # remaining square root gives these two equivalent integer inequalities.
    internal_excess = 2 * ((internal_size + 1) // 2 - 1)
    assert internal_excess >= 0
    remainder = 2 * (n - endpoint_dot) - internal_excess**2 * D
    assert remainder >= 0
    assert remainder**2 >= 8 * internal_excess**2 * D * (n + endpoint_dot)


def main():
    arcs = pairs = domain_rows = internal_rows = 0
    nontrivial_domains = three_row_domains = 0
    for n in range(1, MAX_N + 1):
        points = circle_points(n)
        if len(points) < 2:
            continue
        for block in cyclic_short_arcs(points):
            arcs += 1
            endpoint_dot = block[0][0] * block[-1][0] + block[0][1] * block[-1][1]
            block_set = set(block)
            for first, second in combinations(block, 2):
                numerator = gmul(first, second)
                common = ggcd(numerator, (n, 0))
                reduced_numerator = gdivexact(numerator, common)
                denominator = gdivexact((n, 0), common)
                assert gnorm(ggcd(reduced_numerator, denominator)) == 1
                D = gnorm(denominator)
                assert gnorm(reduced_numerator) == D and n % D == 0

                # The reduced numerator is a unit times conjugate(denominator),
                # so d=conjugate(denominator) is the domain generator.
                unit = gdivexact(reduced_numerator, gconj(denominator))
                assert gnorm(unit) == 1

                domain = []
                internal = []
                for z in block:
                    by_domain = divisible(denominator, gconj(z))
                    reflected_numerator = gmul(reduced_numerator, gconj(z))
                    by_fraction = divisible(denominator, reflected_numerator)
                    assert by_domain == by_fraction
                    if not by_domain:
                        continue
                    reflected = gdivexact(reflected_numerator, denominator)
                    assert gnorm(reflected) == n
                    domain.append(z)
                    if reflected in block_set:
                        internal.append(z)

                assert first in domain and second in domain
                assert first in internal and second in internal

                check_exact_capacity(n, endpoint_dot, D, len(domain), len(internal))

                if D > 1:
                    assert D >= 5
                    nontrivial_domains += 1
                if len(domain) >= 3:
                    three_row_domains += 1

                pairs += 1
                domain_rows += len(domain)
                internal_rows += len(internal)

    assert pairs and nontrivial_domains and three_row_domains
    print(
        f"PASS: {arcs} actual arcs and {pairs} anchor pairs on circles N<={MAX_N}."
    )
    print(
        "PASS: exact reduced domains and reflections; "
        f"{domain_rows} domain-row and {internal_rows} internal-row incidences."
    )
    print(
        "PASS: both stronger exact sagitta inequalities; "
        f"{nontrivial_domains} nonunit and {three_row_domains} three-row domains."
    )


if __name__ == "__main__":
    main()
