"""Exact checks for weighted_rank_one_power_class_packing.md.

Run from the repository root with

    python3 docs/check_weighted_rank_one_power_class_packing.py

The checker uses doubled circle coordinates, so antipodal interval
intersections and all transportation weights remain exact Fractions.
"""

from fractions import Fraction
from itertools import combinations, product
from math import gcd, lcm


def interval_parts(start, length, circumference):
    """Split one half-open circular interval into ordinary intervals."""
    start %= circumference
    end = start + length
    if end <= circumference:
        return ((start, end),)
    return ((start, circumference), (0, end - circumference))


def intersection_length(parts_a, parts_b):
    return sum(
        max(0, min(a1, b1) - max(a0, b0))
        for a0, a1 in parts_a
        for b0, b1 in parts_b
    )


def antipodal_transport(classes):
    """Return y_jh from contiguous class intervals with demands a_j."""
    weights = [weight for cls in classes for weight in cls]
    labels = [label for label, cls in enumerate(classes) for _ in cls]
    total = sum(weights)
    circumference = 2 * total

    intervals = []
    cursor = 0
    for weight in weights:
        intervals.append(interval_parts(cursor, 2 * weight, circumference))
        cursor += 2 * weight
    assert cursor == circumference

    directed = [[Fraction(0) for _ in weights] for _ in weights]
    for j, interval_j in enumerate(intervals):
        shifted_j = tuple(
            piece
            for start, end in interval_j
            for piece in interval_parts(start + total, end - start, circumference)
        )
        for h, interval_h in enumerate(intervals):
            directed[j][h] = Fraction(
                intersection_length(shifted_j, interval_h), 2
            )

    for j in range(len(weights)):
        assert sum(directed[j]) == weights[j]
        for h in range(len(weights)):
            assert directed[j][h] == directed[h][j]
            if labels[j] == labels[h]:
                assert directed[j][h] == 0
    assert sum(
        directed[j][h]
        for j in range(len(weights))
        for h in range(j + 1, len(weights))
    ) == Fraction(total, 2)
    return weights, labels, directed


def audit_case(classes):
    masses = [sum(cls) for cls in classes]
    total = sum(masses)
    assert max(masses) * 2 <= total
    weights, labels, y = antipodal_transport(classes)
    class_gcds = []
    for cls in classes:
        value = 0
        for weight in cls:
            value = gcd(value, weight)
        class_gcds.append(value)

    loads = [Fraction(0) for _ in weights]
    amplified = Fraction(0)
    for j in range(len(weights)):
        for h in range(j + 1, len(weights)):
            if not y[j][h]:
                continue
            aj, ah = weights[j], weights[h]
            common = gcd(aj, ah)
            pair_lcm = lcm(aj, ah)
            certificate_weight = y[j][h] / pair_lcm
            loads[j] += certificate_weight * (ah // common)
            loads[h] += certificate_weight * (aj // common)

            dc = class_gcds[labels[j]]
            dd = class_gcds[labels[h]]
            power = pair_lcm // lcm(dc, dd)
            assert power >= 1
            assert (pair_lcm // dc) % power == 0
            assert (pair_lcm // dd) % power == 0
            amplified += certificate_weight * power

    assert loads == [Fraction(1) for _ in weights]
    expected = sum(
        y[j][h] / lcm(class_gcds[labels[j]], class_gcds[labels[h]])
        for j in range(len(weights))
        for h in range(j + 1, len(weights))
    )
    assert amplified == expected

    h_value = max(
        (class_gcds[c] + class_gcds[d]) // gcd(class_gcds[c], class_gcds[d])
        for c in range(len(classes))
        for d in range(c + 1, len(classes))
    )
    effective_count = sum(
        Fraction(weight, class_gcds[label])
        for weight, label in zip(weights, labels)
    )
    assert amplified >= effective_count / h_value
    return len(weights), amplified


def contiguous_partitions(values):
    """All partitions obtained by cutting an ordered nonempty tuple."""
    if len(values) == 1:
        yield (values,)
        return
    for cuts in product((False, True), repeat=len(values) - 1):
        classes = []
        start = 0
        for index, cut in enumerate(cuts, 1):
            if cut:
                classes.append(values[start:index])
                start = index
        classes.append(values[start:])
        yield tuple(classes)


def check_small_transports():
    cases = 0
    for size in range(2, 7):
        for values in product(range(1, 4), repeat=size):
            for classes in contiguous_partitions(values):
                if len(classes) < 2:
                    continue
                masses = [sum(cls) for cls in classes]
                if 2 * max(masses) > sum(masses):
                    continue
                audit_case(classes)
                cases += 1
    return cases


def check_signed_exponent_roots():
    # Each class gets a distinct prime-coordinate basis vector.  Formula
    # nu_j=(a_j/d_C)t_C realizes exact collision classes.
    examples = (
        ((2, 4), (3, 3, 3, 3)),
        ((5, 10, 15), (6, 9, 12, 18, 30)),
        ((1, 4, 7), (2, 2, 6), (3, 9)),
    )
    pairs = 0
    for classes in examples:
        class_gcds = []
        records = []
        for label, cls in enumerate(classes):
            dc = 0
            for weight in cls:
                dc = gcd(dc, weight)
            class_gcds.append(dc)
            root_vector = [0] * len(classes)
            root_vector[label] = label + 1
            for weight in cls:
                valuation = tuple((weight // dc) * value for value in root_vector)
                records.append((label, weight, valuation))

        for j in range(len(records)):
            cj, aj, nuj = records[j]
            for h in range(j + 1, len(records)):
                ch, ah, nuh = records[h]
                common = gcd(aj, ah)
                exponents = tuple(
                    (ah // common) * x - (aj // common) * z
                    for x, z in zip(nuj, nuh)
                )
                if cj == ch:
                    assert all(value == 0 for value in exponents)
                else:
                    power = lcm(aj, ah) // lcm(class_gcds[cj], class_gcds[ch])
                    assert any(exponents)
                    assert all(value % power == 0 for value in exponents)
                pairs += 1
    return pairs


def check_bounded_height_constant(max_b=30):
    for bound in range(1, max_b + 1):
        actual = max(
            Fraction(c + d, gcd(c, d))
            for c in range(1, bound + 1)
            for d in range(1, bound + 1)
        )
        expected = Fraction(2 if bound == 1 else 2 * bound - 1)
        assert actual == expected


def check_sharp_families(max_b=30):
    for bound in range(2, max_b + 1):
        for scale in (1, 2, 5):
            left_size = scale * (bound - 1)
            right_size = scale * bound
            left_mass = left_size * bound
            right_mass = right_size * (bound - 1)
            assert left_mass == right_mass
            size = left_size + right_size
            amplified = Fraction(left_mass, lcm(bound, bound - 1))
            assert size == scale * (2 * bound - 1)
            assert amplified == scale


def subset_sums(weights):
    values = set()
    for size in range(len(weights) + 1):
        for chosen in combinations(range(len(weights)), size):
            values.add(sum(weights[index] for index in chosen))
    return values


def check_two_class_integer_lines():
    examples = (
        ((2, 4), (3, 3, 3, 3)),
        ((5, 10, 15), (6, 9, 12, 18, 30)),
        ((1, 4, 7), (2, 2, 6)),
    )
    cohorts = 0
    for left, right in examples:
        dc = 0
        dd = 0
        for weight in left:
            dc = gcd(dc, weight)
        for weight in right:
            dd = gcd(dd, weight)
        left_sums = subset_sums(tuple(weight // dc for weight in left))
        right_sums = subset_sums(tuple(weight // dd for weight in right))
        rows_by_load = {}
        for pc, pd in product(left_sums, right_sums):
            rows_by_load.setdefault(dc * pc + dd * pd, []).append((pc, pd))

        common = gcd(dc, dd)
        for rows in rows_by_load.values():
            base_c, base_d = rows[0]
            orbit = []
            for pc, pd in rows:
                assert dc * (pc - base_c) + dd * (pd - base_d) == 0
                assert (pc - base_c) % (dd // common) == 0
                assert (pd - base_d) % (dc // common) == 0
                nc = (pc - base_c) // (dd // common)
                nd = -(pd - base_d) // (dc // common)
                assert nc == nd
                orbit.append(nc)
            if any(orbit):
                delta = 0
                for value in orbit:
                    delta = gcd(delta, abs(value))
                normalized = [value // delta for value in orbit]
                assert max(normalized) - min(normalized) >= len(set(normalized)) - 1
            cohorts += 1
    return cohorts


def gaussian_mul(left, right):
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def gaussian_conj(value):
    return value[0], -value[1]


def gaussian_norm(value):
    return value[0] * value[0] + value[1] * value[1]


def gaussian_power(value, exponent):
    answer = (1, 0)
    while exponent:
        if exponent % 2:
            answer = gaussian_mul(answer, value)
        value = gaussian_mul(value, value)
        exponent //= 2
    return answer


def check_gaussian_orbit_brackets():
    units = ((1, 0), (-1, 0), (0, 1), (0, -1))
    primitive_bases = ((2, 1), (3, 2), (4, 1), (5, 2))
    cases = 0
    for base, exponent, left_unit, right_unit in product(
        primitive_bases, range(1, 8), units, units
    ):
        forward = gaussian_mul(left_unit, gaussian_power(base, exponent))
        backward = gaussian_mul(
            right_unit, gaussian_power(gaussian_conj(base), exponent)
        )
        assert gaussian_norm(forward) == gaussian_norm(backward)
        if forward == backward:
            continue
        difference = forward[0] - backward[0], forward[1] - backward[1]
        squared_distance = gaussian_norm(difference)
        assert squared_distance >= 2 and squared_distance % 2 == 0
        cases += 1
    return cases


def main():
    transport_cases = check_small_transports()
    exponent_pairs = check_signed_exponent_roots()
    check_bounded_height_constant()
    check_sharp_families()
    orbit_cohorts = check_two_class_integer_lines()
    orbit_brackets = check_gaussian_orbit_brackets()
    print(
        "verified",
        transport_cases,
        "exact antipodal transports and certificate load systems",
    )
    print("verified", exponent_pairs, "same-class collisions/cross-class power divisibilities")
    print("verified H_B through B=30 and the sharp two-class families at three scales")
    print("verified", orbit_cohorts, "two-class weighted subset-sum orbit cohorts")
    print("verified", orbit_brackets, "noncoincident equal-norm Gaussian orbit brackets")
    print("Roth separation and the whole-row gcd class cap are prose inputs, not numerical tests")


if __name__ == "__main__":
    main()
