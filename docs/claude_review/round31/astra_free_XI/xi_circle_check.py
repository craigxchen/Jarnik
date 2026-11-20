"""Exact minimum Euclidean chord on the squarefree split parent in [P, 2P]."""

import argparse
from functools import cmp_to_key
from math import isqrt, prod


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, isqrt(n) + 1))


def sum_of_two_squares(p):
    for x in range(1, isqrt(p) + 1):
        y = isqrt(p - x * x)
        if x * x + y * y == p:
            return x, y
    raise ArithmeticError(f"missing split factor for {p}")


def first_quadrant(x, y):
    while x <= 0 or y < 0:
        x, y = -y, x
    return x, y


def compare_directions(a, b):
    return a[1] * b[0] - a[0] * b[1]


def squared_distance(a, b):
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2


def check_parent(P):
    primes = [p for p in range(P, 2 * P + 1) if p % 4 == 1 and is_prime(p)]
    if not primes:
        raise ValueError("the interval must contain at least one split prime")
    q = prod(primes)
    points = [(1, 0)]
    for p in primes:
        a, b = sum_of_two_squares(p)
        new_points = []
        for x, y in points:
            new_points.append((a * x - b * y, b * x + a * y))
            new_points.append((a * x + b * y, a * y - b * x))
        points = new_points
    del new_points

    # Every unit orbit occurs once by unique factorization in Z[i].
    for j, (x, y) in enumerate(points):
        if x * x + y * y != q:
            raise ArithmeticError("incorrect norm")
        points[j] = first_quadrant(x, y)
    points.sort(key=cmp_to_key(compare_directions))

    # The four quadrants are quarter-turn copies. Include the last gap
    # of this quadrant, running from its last point to i times its first.
    first = points[0]
    next_quadrant_first = (-first[1], first[0])
    endpoints = (points[-1], next_quadrant_first)
    minimum = squared_distance(*endpoints)
    previous = first
    for j in range(1, len(points)):
        current = points[j]
        if compare_directions(previous, current) >= 0:
            raise ArithmeticError("duplicate orbit or incorrect ordering")
        distance = squared_distance(previous, current)
        if distance < minimum:
            minimum = distance
            endpoints = previous, current
        previous = current

    floor_ratio = isqrt(isqrt(minimum ** 2 // q))
    print(f"P={P} s={len(primes)} orbits={len(points)} q={q}")
    print(f"min_chord_squared={minimum} endpoints={endpoints}")
    print(f"floor(min_chord/q^(1/4))={floor_ratio}")
    print(f"min_chord>8*q^(1/4): {minimum ** 2 > 8 ** 4 * q}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("P", type=int, nargs="*", default=[200, 250])
    for lower in parser.parse_args().P:
        check_parent(lower)
