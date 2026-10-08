"""Exact finite search for resonant cyclotomic prefix fibres.

For odd e>1, a layer of width W_e chooses a_e in [0,W_e] between
Phi_e(z) and Phi_e(-z)=Phi_{2e}(z).  The base layer has even width W_1
and even choice a_1.  The total reciprocal degree is

    L = W_1 + sum_{odd e>1} phi(e) W_e.

This checker exhausts all profiles with L<=40.  It groups choices by
their exact coefficient prefix of orders <ceil(L/4), using Ramanujan
moments rather than expanding the full polynomials.  It also imposes the
necessary weighted-distance condition d>=2ceil(L/4) on every pair in
every fibre; no rows are discarded.

The search is finite evidence only; it is not a uniform theorem.
"""

from collections import defaultdict
from itertools import combinations
from math import gcd


CAP = 40


def phi(n: int) -> int:
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            result = result // p * (p - 1)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        result = result // n * (n - 1)
    return result


def mobius(n: int) -> int:
    result = 1
    p = 2
    while p * p <= n:
        if n % (p * p) == 0:
            return 0
        if n % p == 0:
            result = -result
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        result = -result
    return result


def ramanujan(e: int, n: int) -> int:
    """c_e(n)=mu(e/g)*phi(e)/phi(e/g), g=gcd(e,n)."""
    g = gcd(e, n)
    d = e // g
    mu = mobius(d)
    return 0 if mu == 0 else mu * phi(e) // phi(d)


def profile_fibres(profile):
    """Return exact moment-signature fibres for one W profile."""
    degree = sum(phi(e) * width for e, width in profile.items())
    cutoff = (degree + 3) // 4
    orders = tuple(range(1, cutoff, 2))
    factors = tuple(sorted(profile))
    weights = tuple(phi(e) for e in factors)
    choices = tuple(
        tuple(range(0, profile[e] + 1, 2)) if e == 1
        else tuple(range(profile[e] + 1))
        for e in factors
    )
    columns = tuple(
        tuple(ramanujan(e, n) for n in orders) for e in factors
    )
    groups = defaultdict(list)
    row = [0] * len(factors)

    def visit(j, signature):
        if j == len(factors):
            groups[signature].append(tuple(row))
            return
        column = columns[j]
        for exponent in choices[j]:
            row[j] = exponent
            visit(j + 1, tuple(x + exponent * y
                                for x, y in zip(signature, column)))

    visit(0, (0,) * len(orders))
    return degree, cutoff, orders, factors, weights, groups


def multiply_truncated(left, right, length):
    result = [0] * length
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            if i + j < length:
                result[i + j] += x * y
    return result


def check_four_row_regression() -> None:
    """Astra's L=12 four-row trade, with a nonzero common prefix."""
    factors = {
        3: ([1, 1, 1], [1, -1, 1]),
        5: ([1, 1, 1, 1, 1], [1, -1, 1, -1, 1]),
        9: ([1, 0, 0, 1, 0, 0, 1], [1, 0, 0, -1, 0, 0, 1]),
    }
    rows = ((1, 0, 0), (1, 0, 1), (0, 1, 0), (0, 1, 1))
    prefixes = []
    for exponents in rows:
        polynomial = [1]
        for e, orientation in zip((3, 5, 9), exponents):
            polynomial = multiply_truncated(
                polynomial, factors[e][0 if orientation else 1], 3
            )
        prefixes.append(tuple(polynomial))
    assert prefixes == [(1, 0, 1)] * 4
    weights = (phi(3), phi(5), phi(9))
    distances = [
        sum(w * abs(x - y) for w, x, y in zip(weights, left, right))
        for left, right in combinations(rows, 2)
    ]
    assert min(distances) == 6
    assert 2 * ((12 + 3) // 4) == 6
    assert all(
        sum(a * ramanujan(e, 1) for e, a in zip((3, 5, 9), row)) == -1
        for row in rows
    )

    # A nonzero common odd coefficient represents a moving arc center.
    rows = ((0, 0, 0, 0), (2, 0, 1, 1),
            (2, 1, 0, 1), (2, 1, 1, 0))
    prefixes = []
    for base, *prime_counts in rows:
        polynomial = [1]
        for _ in range(base):
            polynomial = multiply_truncated(polynomial, [1, -1], 10)
        for _ in range(2 - base):
            polynomial = multiply_truncated(polynomial, [1, 1], 10)
        for p, orientation in zip((11, 13, 17), prime_counts):
            factor = [1 if orientation else (-1) ** j for j in range(p)]
            polynomial = multiply_truncated(polynomial, factor, 10)
        prefixes.append(tuple(polynomial))
    assert len(set(prefixes)) == 1
    assert prefixes[0][1] == -1


def exhaustive_profiles():
    """Yield all W profiles with weighted degree <=CAP."""
    # For odd e=prod p^a, phi(e)^2/e is multiplicative with local
    # factors p^(a-2)(p-1)^2 >= 1. Hence phi(e)^2>=e, so
    # phi(e)<=CAP implies e<=CAP^2.
    layers = []
    for e in range(3, CAP * CAP + 1, 2):
        weight = phi(e)
        if weight <= CAP:
            assert weight * weight >= e
            layers.append((e, weight))

    profile = {}

    def visit(index, remaining):
        if index == len(layers):
            for base_width in range(0, remaining + 1, 2):
                result = dict(profile)
                if base_width:
                    result[1] = base_width
                degree = sum(phi(e) * w for e, w in result.items())
                if degree:
                    yield result
            return
        e, weight = layers[index]
        for width in range(remaining // weight + 1):
            if width:
                profile[e] = width
            yield from visit(index + 1, remaining - width * weight)
        profile.pop(e, None)

    yield from visit(0, CAP)


def run_search() -> None:
    profile_count = 0
    orientation_count = 0
    signature_count = 0
    checked_pairs = 0
    max_fibre = 0
    max_fibre_fixture = None
    max_by_degree = {}

    for profile in exhaustive_profiles():
        degree, cutoff, orders, factors, weights, groups = profile_fibres(profile)
        profile_count += 1
        orientation_count += sum(len(rows) for rows in groups.values())
        signature_count += len(groups)
        local_max = max(len(rows) for rows in groups.values())
        if local_max > max_fibre:
            max_fibre = local_max
            sig, rows = next((sig, rows) for sig, rows in groups.items()
                             if len(rows) == local_max)
            max_fibre_fixture = (degree, profile, orders, sig, rows)
        max_by_degree[degree] = max(max_by_degree.get(degree, 0), local_max)

        threshold = 2 * cutoff
        for rows in groups.values():
            high_coordinates = [j for j, e in enumerate(factors) if e >= cutoff]
            projections = {tuple(row[j] for j in high_coordinates) for row in rows}
            assert len(projections) == len(rows)
            for left, right in combinations(rows, 2):
                checked_pairs += 1
                distance = sum(
                    weight * abs(x - y)
                    for weight, x, y in zip(weights, left, right)
                )
                assert distance >= threshold, (
                    degree, profile, orders, left, right, distance, threshold
                )

    assert profile_count == 20886
    assert orientation_count == 1506841
    assert signature_count == 1437928
    assert max_fibre == 4
    print(f"degree cap: {CAP}")
    print(f"active odd layer indices: {sum(1 for e in range(3, CAP*CAP+1, 2) if phi(e)<=CAP)}")
    print(f"profiles: {profile_count}")
    print(f"orientation vectors: {orientation_count}")
    print(f"distinct signatures across profiles: {signature_count}")
    print(f"within-fibre pairs checked: {checked_pairs}")
    print(f"maximum raw fibre: {max_fibre}; fixture: {max_fibre_fixture}")
    print("max raw fibre by L:", sorted(max_by_degree.items()))


if __name__ == "__main__":
    check_four_row_regression()
    run_search()
