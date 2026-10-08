"""Independent projective-root and contact audit; exact rational arithmetic."""
from fractions import Fraction as F
from itertools import combinations, product

PARAMETERS = (F(0), F(-240, 5329), F(8, 257))
ROOTS = (F(-1), F(3), F(4), F(5))
TARGETS = (F(-1, 2), F(9, 2), F(16, 3), F(25, 4))
LABELS = frozenset(("0", "1", "inf", "a", "b", "c"))


def multiply(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def polynomials(lam):
    return (-60 * lam, -73 * lam, F(1)), (-1 - 38 * lam, 1 - 11 * lam, lam)


def evaluate(p, x):
    return sum(c * x**i for i, c in enumerate(p))


def normalized(x):
    return None if x == 4 else (x + 1) / (16 - 4 * x)


def boundary(side):
    side = frozenset(side)
    other = LABELS - side
    assert min(len(side), len(other)) >= 2
    return min((tuple(sorted(side)), tuple(sorted(other))), key=lambda x: (len(x), x))


def main():
    determinant = [F(1)]
    for root in ROOTS:
        determinant = multiply(determinant, [-root, F(1)])
    for lam, mu in combinations(PARAMETERS, 2):
        p, q = polynomials(lam)
        r, s = polynomials(mu)
        assert [x - y for x, y in zip(multiply(p, s), multiply(r, q))] == [
            (mu - lam) * x for x in determinant
        ]
    # At infinity the numerator leading coefficient is 1, and the
    # homogeneous pair cross-products have nonzero leading coefficient.
    assert len({None if lam == 0 else 1 / lam for lam in PARAMETERS}) == 3
    for lam in PARAMETERS:
        p, q = polynomials(lam)
        # Any finite common factor must vanish at a determinant root.
        assert all(evaluate(q, x) != 0 for x in ROOTS)
        assert p[2] == 1
        # The derivative numerator has degree two, hence no critical
        # point at infinity; there is no double finite pole either.
        assert p[2] * q[1] - p[1] * q[2] != 0
        assert q[2] == 0 or q[1]**2 - 4 * q[2] * q[0] != 0

    points, contacts, clusters = set(), {}, []
    for root, target in zip(ROOTS, TARGETS):
        common = normalized(root)
        alternatives, slopes = [], []
        for lam in PARAMETERS:
            p, q = polynomials(lam)
            assert evaluate(p, root) == target * evaluate(q, root)
            a, b, c = (p[i] - target * q[i] for i in (2, 1, 0))
            assert a != 0  # Both roots are finite in the original x-chart.
            alt = -b / a - root
            assert a * alt**2 + b * alt + c == 0 and alt != root
            alternatives.append(normalized(alt))
            d = ((p[1] + 2 * p[2] * root) * evaluate(q, root)
                 - evaluate(p, root) * (q[1] + 2 * q[2] * root)) / evaluate(q, root)**2
            assert d != 0
            chart_derivative = (-20 / (root + 1)**2 if root == 4
                                else 20 / (16 - 4 * root)**2)
            slopes.append(chart_derivative / d)
        assert len(set(alternatives)) == 3
        assert not set(alternatives) & {F(0), F(1), None, common}
        assert len(set(slopes)) == 3 and all(slopes)
        anchor = {F(0): "0", F(1): "1", None: "inf"}.get(common)
        fiber_contacts = 0
        for choices in product((False, True), repeat=3):
            coords = tuple(common if bit else alt for bit, alt in zip(choices, alternatives))
            assert coords not in points
            points.add(coords)
            selected = {label for label, bit in zip(("a", "b", "c"), choices) if bit}
            if anchor is not None:
                selected.add(anchor)
            if len(selected) >= 2:
                key = boundary(selected)
                assert key not in contacts
                contacts[key] = coords
                clusters.append(frozenset(selected))
                fiber_contacts += 1
        assert fiber_contacts == (7 if anchor is not None else 4)
    all_boundaries = {boundary(side) for n in (2, 3, 4) for side in combinations(sorted(LABELS), n)}
    assert len(points) == 32
    assert set(contacts) == all_boundaries and len(contacts) == 25
    assert len(set(contacts.values())) == 25
    for left, right in combinations(sorted(LABELS), 2):
        moving_count = sum(label in ("a", "b", "c") for label in (left, right))
        assert sum({left, right} <= cluster for cluster in clusters) == 4 * moving_count
    print("Exact projective determinant, degree-two maps, and no extra pole contacts verified.")
    print("32 distinct rational points in four fibers; 25 distinct boundary points verified.")
    print("All common-root inverse slopes are nonzero and pairwise distinct.")
    print("Pair determinant divisor degrees are 0, 4, and 8 for zero, one, and two moving labels.")


if __name__ == "__main__":
    main()
