"""Exact certificate for quartic_extension_maximality.md; no dependencies."""

from fractions import Fraction as Q
from itertools import combinations


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(p, q):
    return trim([
        (p[j] if j < len(p) else Q(0))
        + (q[j] if j < len(q) else Q(0))
        for j in range(max(len(p), len(q)))
    ])


def mul(p, q):
    result = [Q(0)] * (len(p) + len(q) - 1)
    for j, a in enumerate(p):
        for k, b in enumerate(q):
            result[j + k] += a * b
    return trim(result)


def remainder(p, q):
    p, q = trim(p), trim(q)
    assert q != [0]
    while p != [0] and len(p) >= len(q):
        shift, coefficient = len(p) - len(q), p[-1] / q[-1]
        p = add(p, [Q(0)] * shift + [-coefficient * a for a in q])
    return p


def vector(numerators, denominator=1):
    p = tuple(Q(n, denominator) for n in numerators)
    return p + (Q(0),) * (9 - len(p))


ANCHORS = [
    vector([112, -152, 121, -48, 9], 104),
    vector([26, -421, 8, 21, -18], 442),
    vector([189, -944, 162, 144, -27], 1088),
]
EXTRA = [
    vector([-5744, 4024, -2277, 576, -108], 952),
    vector([-13966, 11661, -5328, 1539, -162], 3978),
    vector([37193, -4728, -531, 1728, -324], 19656),
]
ROOTS = [
    [(0, 1), (2, 2), (Q(8, 3), -1), (Q(2, 3), -2)],
    [(0, 1), (2, 2), (-Q(8, 3), -Q(1, 3)), (Q(11, 6), -Q(8, 3))],
    [(0, 1), (Q(8, 3), -1), (-Q(8, 3), -Q(1, 3)), (Q(16, 3), Q(1, 3))],
]


def directions(anchor, roots):
    quadratics = [
        [Q(a) ** 2 + Q(b) ** 2, -2 * Q(a), Q(1)]
        for a, b in roots
    ]
    assert len({tuple(p) for p in quadratics}) == 4
    assert all(b != 0 for _, b in roots)
    products = []
    for size in range(5):
        for chosen in combinations(quadratics, size):
            p = [Q(1)]
            for factor in chosen:
                p = mul(p, factor)
            products.append(vector(p))
    norm = add(mul(trim(anchor), trim(anchor)), [Q(1)])
    assert norm == [anchor[4] ** 2 * a for a in products[-1]]
    assert len(set(products)) == 16
    return products


def intersection(a, p, b, q):
    """Return a unique line intersection, or None for parallel/disjoint lines."""
    delta = tuple(bb - aa for aa, bb in zip(a, b))
    for j, k in combinations(range(9), 2):
        determinant = -p[j] * q[k] + p[k] * q[j]
        if determinant:
            lam = (-delta[j] * q[k] + delta[k] * q[j]) / determinant
            mu = (p[j] * delta[k] - p[k] * delta[j]) / determinant
            if all(lam * p[l] - mu * q[l] == delta[l] for l in range(9)):
                return tuple(a[l] + lam * p[l] for l in range(9))
            return None
    return None


def compatible(candidate, anchor):
    if candidate == anchor:
        return True
    difference = [a - b for a, b in zip(candidate, anchor)]
    numerator = add(mul(trim(candidate), trim(anchor)), [Q(1)])
    return remainder(numerator, difference) == [0]


def main():
    u = [a - b for a, b in zip(ANCHORS[1], ANCHORS[0])]
    v = [a - b for a, b in zip(ANCHORS[2], ANCHORS[0])]
    assert any(u[j] * v[k] != u[k] * v[j] for j, k in combinations(range(9), 2))
    all_directions = [directions(a, r) for a, r in zip(ANCHORS, ROOTS)]
    candidates = set()
    for a, b in combinations(range(3), 2):
        for p in all_directions[a]:
            for q in all_directions[b]:
                candidate = intersection(ANCHORS[a], p, ANCHORS[b], q)
                if candidate is not None:
                    candidates.add(candidate)
    assert candidates == set(ANCHORS + EXTRA)
    assert [tuple(compatible(f, a) for a in ANCHORS) for f in EXTRA] == [
        (True, True, False), (True, False, True), (False, True, True)
    ]
    survivors = {f for f in candidates if all(compatible(f, a) for a in ANCHORS)}
    assert survivors == set(ANCHORS)
    print("Passed: 768 exact line-pair checks; six candidates; only three original anchors survive.")


if __name__ == "__main__":
    main()
