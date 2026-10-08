"""Exact rational checks of actual-point inversions and their least-radius cost."""

from fractions import Fraction as Q
from itertools import combinations
from math import atan, gcd, log


def chart(x, length):
    denominator = x * x + length * length
    return ((x * x - length * length) / denominator,
            2 * x * length / denominator)


def invert(point, center, scale):
    dx, dy = point[0] - center[0], point[1] - center[1]
    square = dx * dx + dy * dy
    return (center[0] + scale * dx / square,
            center[1] + scale * dy / square)


def edge_norm(x, length):
    a, b = x.numerator, length * x.denominator
    common = gcd(a, b)
    a, b = a // common, b // common
    epsilon = 2 if a % 2 and b % 2 else 1
    return (a * a + b * b) // epsilon


def least_square(coordinates, length):
    edges = list(map(Q, coordinates))
    edges.extend(Q(x * y + length * length, x - y)
                 for x, y in combinations(coordinates, 2))
    value = 1
    for x in edges:
        contribution = edge_norm(x, length)
        value = value // gcd(value, contribution) * contribution
    return value


def pell(index):
    u, v = 1, 0
    for _ in range(index):
        u, v = 9 * u + 20 * v, 4 * u + 9 * v
    assert u * u - 5 * v * v == 1
    return u, v


def main():
    identities = 0
    for length in (Q(1), Q(3), Q(24)):
        for x in (Q(-7), Q(-1, 3), Q(2, 5), Q(1), Q(13)):
            for t in (Q(1, 2), Q(1), Q(2), Q(3), Q(7, 2)):
                source = chart(x, length)
                first = invert(source, (Q(1), Q(0)), t)
                assert first == (1 - t / 2, t * x / (2 * length))
                second = invert(first, (Q(-1), Q(0)), 4 - t)
                k = (4 - t) / t
                assert second == chart(k * length * length / x, length)
                reanchored = (-second[0], second[1])
                assert reanchored == chart(x / k, length)
                assert sum(z * z for z in second) == 1
                identities += 1

    # The four-point selective-map exception to the five-point rigidity count.
    for x, expected in ((Q(1), Q(1)), (Q(-1), Q(-1)),
                        (Q(2), Q(-2)), (Q(1, 2), Q(-1, 2))):
        assert (5 * x - 4) / (5 - 4 * x) == expected

    observations = []
    for index in (1, 11, 21):
        u, v = pell(index)
        xs = (90 * u * v + 30 * v * v + 3,
              96 * u * v, 160 * v * v + 128 * u * v + 16)
        a = 1 + 10 * v * v + 2 * u * v
        b = 1 + 10 * v * v - 2 * u * v
        c = 13 + 130 * v * v + 38 * u * v
        source = least_square(xs, 24)
        target = least_square(tuple(3 * x for x in xs), 24)
        assert source == 5 * a * b * c
        observations.append((len(str(source)), len(str(target))))
        if index == 1:
            assert xs == (3723, 3456, 7184)
            assert source == 358853785
            assert target == 2085984659911668125
            source_span = 2 * source ** 0.25 * atan(24 / min(xs))
            target_span = 2 * target ** 0.25 * atan(8 / min(xs))
            assert 1.911 < source_span < 1.913
            assert 175.943 < target_span < 175.945
        else:
            ratio = log(target) / log(source)
            assert 2 < ratio < 2.015
    assert observations == [(9, 19), (84, 169), (160, 320)]
    print("PASS: %d exact inversion identities; four-point projective fixture; "
          "three exact Pell all-edge radius computations." % identities)


if __name__ == "__main__":
    main()
