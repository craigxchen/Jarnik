"""Exact Fourier-gap audit and the full M=32 character-height minimum."""

from itertools import product

from check_walsh_nonlinear_coset_height_obstruction import (
    B, M, LABELS, exhaustive_counts, physical_l1,
)


def parity(n):
    return bin(n).count("1") & 1


def signed_coset_character(values):
    support = [x for x, value in enumerate(values) if value]
    if not support:
        return False
    origin = support[0]
    translated = {x ^ origin for x in support}
    if any(x ^ y not in translated for x in translated for y in translated):
        return False
    normalized = {
        x: values[x ^ origin] * values[origin] for x in translated
    }
    return all(normalized[x ^ y] == normalized[x] * normalized[y]
               for x in translated for y in translated)


def raw_fourier_l1(values):
    return sum(abs(sum(value * (-1 if parity(a & x) else 1)
                       for x, value in enumerate(values)))
               for a in range(len(values)))


def check_fourier_gap():
    n = 8
    counts = {"zero": 0, "coset": 0, "gap": 0}
    for values in product((-1, 0, 1), repeat=n):
        if not any(values):
            counts["zero"] += 1
            continue
        total = raw_fourier_l1(values)
        if signed_coset_character(values):
            assert total == n
            counts["coset"] += 1
        else:
            assert 2 * total >= 3 * n
            counts["gap"] += 1
    # The three-point indicator attains the noncoset gap.
    assert raw_fourier_l1((1, 1, 1, 0, 0, 0, 0, 0)) == 12
    assert sum(counts.values()) == 3 ** n
    return counts


def physical_integer_height(c, b, labels):
    """Exact copies: reserve one separate copy for each assigned flip."""
    m = len(c)
    result = 0
    for a in range(1, m):
        twice = sum(cx * (-1 if parity(a & x) else 1)
                    for x, cx in enumerate(c))
        assert twice % 2 == 0
        baseline = twice // 2
        assigned = [x for x in range(m) if labels[x] == a]
        assert len(assigned) <= b
        result += (b - len(assigned)) * abs(baseline)
        result += sum(abs(baseline - c[x] * (-1 if parity(a & x) else 1))
                      for x in assigned)
    return result


def check_integer_bounds():
    # All zero-sum integer coefficient vectors in a small exact fixture,
    # including coefficients beyond {0,+1,-1}.
    labels = (1, 1, 2, 3)
    count = 0
    for c in product(range(-3, 4), repeat=4):
        if sum(c) or not any(c):
            continue
        amplitude = max(map(abs, c))
        for b in (4, 5, 7):
            height = physical_integer_height(c, b, labels)
            raw_norm = raw_fourier_l1(c)
            assert 2 * height >= b * raw_norm - 2 * sum(map(abs, c))
            if amplitude >= 2 or not signed_coset_character(c):
                assert height >= b * 4 // 2
            count += 1
    return count


def check_full_m32_minimum():
    assert sorted(LABELS) == sorted(tuple(range(1, M)) + (3,))
    counts = exhaustive_counts()
    maximum_matches = [max(counts[d][1]) for d in range(1, 6)]
    assert maximum_matches == [2, 3, 4, 3, 2]
    coset_minima = [B * M // 2 + (1 << d) - 2 * maximum_matches[d - 1]
                   for d in range(1, 6)]
    assert coset_minima == [78, 78, 80, 90, 108]
    assert (3 * B - 4) * M >= 2 * B * M
    assert (B - 2) * M >= B * M // 2
    assert min(coset_minima + [B * M // 2]) == 78
    assert 2 * 78 > B * (M - 1)
    # Direct physical-column construction attains 78.
    assert physical_l1((0, 5, 10, 15), 0, 28) == 78
    return coset_minima


def main():
    counts = check_fourier_gap()
    integer_count = check_integer_bounds()
    minima = check_full_m32_minimum()
    print("PASS: all 6561 F2^3 ternary functions %s; %d integer fixtures; "
          "all M32 coset minima %s; full saturated minimum 78 > 155/2."
          % (counts, integer_count, minima))


if __name__ == "__main__":
    main()
