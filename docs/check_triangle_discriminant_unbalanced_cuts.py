"""Exact local squareclasses for six-point singleton and pair cuts."""

from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd, lcm


def product(values):
    return reduce(lambda x, y: x * y, values, 1)


def central_weights(nodes):
    raw = [Fraction(r * r,
                    product(r - s for j, s in enumerate(nodes) if j != i))
           for i, r in enumerate(nodes)]
    denominator = lcm(*(x.denominator for x in raw))
    values = [int(x * denominator) for x in raw]
    common = reduce(gcd, map(abs, values))
    values = [x // common for x in values]
    assert all(sum(Fraction(u, r ** (-j)) if j < 0 else u * r ** j
                   for u, r in zip(values, nodes)) == 0
               for j in range(-2, 3))
    return values


def valuation(value, prime):
    assert value
    value = abs(value)
    exponent = 0
    while value % prime == 0:
        exponent += 1
        value //= prime
    return exponent


def local_class(value, prime):
    exponent = valuation(value, prime)
    unit = (value // prime ** exponent) % prime
    symbol = pow(unit, (prime - 1) // 2, prime)
    assert symbol in (1, prime - 1)
    return exponent % 2, 1 if symbol == 1 else -1


def discriminant_class(weights, triple):
    subsum = sum(weights[i] for i in triple)
    assert subsum
    return -subsum * product(weights[i] for i in triple)


def check_pair_cut(prime, height):
    nodes = [prime ** height, 2 * prime ** height, 1, 2, 3, 4]
    weights = central_weights(nodes)
    assert [valuation(u, prime) for u in weights] == [height, height, 0, 0, 0, 0]

    # Check the residue formulas (5), including their common unit c.
    ys = nodes[2:]
    poly_derivative = lambda y: product(y - z for z in ys if z != y)
    common_units = {(weights[j] * poly_derivative(nodes[j])) % prime
                    for j in range(2, 6)}
    assert len(common_units) == 1
    common = common_units.pop()
    y_product = product(ys) % prime
    x1, x2 = 1, 2
    eta1 = (weights[0] // prime ** height) % prime
    eta2 = (weights[1] // prime ** height) % prime
    denominator = (x1 - x2) * y_product % prime
    assert eta1 == common * x1 * x1 * pow(denominator, -1, prime) % prime
    assert eta2 == -common * x2 * x2 * pow(denominator, -1, prime) % prime

    triples = list(combinations(range(6), 3))
    for triple in triples:
        crossing_count = len(set(triple) & {0, 1})
        value = discriminant_class(weights, triple)
        if crossing_count in (0, 2):
            assert local_class(value, prime) == (0, 1)
        complement = tuple(i for i in range(6) if i not in triple)
        assert local_class(value, prime) == local_class(
            discriminant_class(weights, complement), prime)

    crossing = [local_class(discriminant_class(weights, triple), prime)
                for triple in triples if len(set(triple) & {0, 1}) == 1]
    # Check the explicit unit squareclass (13) for triples using minority 0.
    for majority_pair in combinations(range(2, 6), 2):
        triple = (0,) + majority_pair
        subsum = sum(weights[i] for i in triple)
        if valuation(subsum, prime):
            continue
        value = discriminant_class(weights, triple)
        assert valuation(value, prime) == height
        complement_pair = [i for i in range(2, 6) if i not in majority_pair]
        difference = (sum(nodes[i] for i in majority_pair)
                      - sum(nodes[i] for i in complement_pair)) % prime
        predicted = -difference * pow(denominator, -1, prime) % prime
        actual = value // prime ** height % prime
        quotient = actual * pow(predicted, -1, prime) % prime
        assert pow(quotient, (prime - 1) // 2, prime) == 1
    if height == 1:
        # All three possibilities really occur; no crossing splitting rule exists.
        assert (1, 1) in crossing or (1, -1) in crossing
        assert (0, 1) in crossing
        assert (0, -1) in crossing
    return weights


def check_singleton_cut(prime):
    nodes = [prime, 1, 2, 3, 4, 7]
    weights = central_weights(nodes)
    assert [valuation(u, prime) for u in weights] == [2, 0, 0, 0, 0, 0]

    triples = list(combinations(range(6), 3))
    classes = []
    for triple in triples:
        value = discriminant_class(weights, triple)
        complement = tuple(i for i in range(6) if i not in triple)
        assert local_class(value, prime) == local_class(
            discriminant_class(weights, complement), prime)
        classes.append(local_class(value, prime))
    assert (0, 1) in classes       # split
    assert (0, -1) in classes      # unramified nonsplit
    assert any(parity for parity, _ in classes)  # ramified
    return weights


def main():
    pair_one = check_pair_cut(13, 1)
    pair_two = check_pair_cut(13, 2)
    singleton = check_singleton_cut(13)
    print("PASS: exact pair heights 1 and 2; four noncrossing local squares")
    print("PASS: crossed pair and singleton fields realize split, nonsplit, and ramified classes")
    print("Witness vectors:", pair_one, pair_two, singleton)


if __name__ == "__main__":
    main()
