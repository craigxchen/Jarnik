"""Bounded exact search for positive six-label cotangent cliques.

The search parametrizes a clique by rational ratios t_i=X_i/L.  For any
five positive integers t_i, clearing the denominators of all ten pair
quotients gives an integral clique.  The test is deliberately finite: it
records examples without making an asymptotic claim.
"""

from itertools import combinations
from math import gcd, lcm

from check_integer_cotangent_normalization import check_clique, primitive_tuple
from check_integer_cotangent_spectral import (
    ONE, div, mul, neg, power, primitive_vector,
)


def scaled_clique(ts):
    """Return the least common denominator scale for the ratios ts."""
    scale = 1
    for x, y in combinations(ts, 2):
        numerator = x * y + 1
        denominator = abs(x - y)
        scale = lcm(scale, denominator // gcd(denominator, numerator))
    return tuple(scale * t for t in ts), scale


def inspect(ts):
    xs, scale = scaled_clique(ts)
    check_clique(xs, scale)

    nodes = [ONE]
    nodes.extend(div((x, -scale), (x, scale)) for x in xs)
    weights = []
    for i, node in enumerate(nodes):
        denominator = ONE
        for j, other in enumerate(nodes):
            if i != j:
                denominator = mul(denominator, add(node, neg(other)))
        weights.append(div(ONE, denominator))

    raw = [mul(weights[i], power(nodes[i], 2)) for i in range(6)]
    central, _ = primitive_vector(raw)
    assert all(z[1] == 0 for z in central)
    central = tuple(z[0] for z in central)

    rows = primitive_tuple(xs, scale)
    products = {}
    collisions = []
    units = ((1, 0), (-1, 0), (0, 1), (0, -1))
    for i, j in combinations(range(6), 2):
        product = mul(rows[i], rows[j])
        if any(mul(unit, product) in products for unit in units):
            collisions.append(((i, j),))
        for unit in units:
            products[mul(unit, product)] = (i, j)

    radius_squared = rows[0][0] ** 2 + rows[0][1] ** 2
    pair_quotients = [
        (x * y + scale * scale) // (x - y)
        for x, y in combinations(xs, 2)
    ]
    return {
        "ts": tuple(ts),
        "xs": xs,
        "scale": scale,
        "radius_squared": radius_squared,
        "min_ratio": xs[0] // scale,
        "max_cotangent": max(max(xs), max(abs(q) for q in pair_quotients)),
        "central": central,
        "collisions": tuple(collisions),
    }


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def main():
    found = []
    tested = 0
    for ts in combinations(range(1, 13), 5):
        data = inspect(ts)
        tested += 1
        if not data["collisions"]:
            found.append(data)

    found.sort(key=lambda d: (d["scale"], d["radius_squared"]))
    print(f"tested {tested} five-tuples; associate-collision-free: {len(found)}")
    for data in found[:3]:
        print(data)


if __name__ == "__main__":
    main()
