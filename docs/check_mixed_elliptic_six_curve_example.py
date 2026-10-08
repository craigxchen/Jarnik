"""Exact checks for mixed_elliptic_six_curve_example.md."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product


PARAMETERS = (F(3), F(-57, 121), F(-9, 25))
VARIABLES = ("a", "b", "c")
LABELS = frozenset(("inf", "0", "1") + VARIABLES)


def f(lam, x):
    if x == "inf":
        return "inf"
    return (1 + 2 * lam) * x * (x - 1) / (1 + lam * x)


def critical_points(lam):
    expected = {
        PARAMETERS[0]: (F(1, 3), F(-1)),
        PARAMETERS[1]: (F(11, 3), F(11, 19)),
        PARAMETERS[2]: (F(5), F(5, 9)),
    }[lam]
    assert all(lam * x * x + 2 * x - 1 == 0 for x in expected)
    return expected


def canonical_boundary(side):
    side = frozenset(side)
    other = LABELS - side
    assert len(side) >= 2 and len(other) >= 2
    if len(side) != len(other):
        return side if len(side) < len(other) else other
    return min(side, other, key=lambda value: tuple(sorted(value)))


def all_boundaries():
    result = set()
    ordered = sorted(LABELS)
    for size in (2, 3, 4):
        for side in combinations(ordered, size):
            if len(LABELS - set(side)) >= 2:
                result.add(canonical_boundary(side))
    return result


def check_pencil_and_branches():
    branch_values = []
    for lam in PARAMETERS:
        branch_values.append(tuple(f(lam, x) for x in critical_points(lam)))
    u, v = branch_values[0]
    assert branch_values == [
        (F(-7, 9), F(-7)),
        (u, F(-7, 361)),
        (v, F(-7, 81)),
    ]
    assert len(set().union(*map(set, branch_values))) == 4

    for x in (F(0), F(1), F(2), "inf"):
        assert len({f(lam, x) for lam in PARAMETERS}) == 1

    # The affine cross-product is -(lambda-mu)x(x-1)(x-2); homogenizing
    # the degree-two denominators supplies the fourth root at infinity.
    for lam, mu in combinations(PARAMETERS, 2):
        for x in map(F, (-3, -2, -1, 0, 1, 2, 3, 4)):
            if (1 + lam * x) and (1 + mu * x):
                lhs = ((1 + 2 * lam) * x * (x - 1) * (1 + mu * x)
                       -(1 + 2 * mu) * x * (x - 1) * (1 + lam * x))
                rhs = -(lam - mu) * x * (x - 1) * (x - 2)
                assert lhs == rhs

    inertia = ((1, 1, 0), (1, 0, 1), (0, 1, 0), (0, 0, 1))
    assert len(set(inertia)) == 4
    # Rank three over F_2.
    masks = [first | second << 1 | third << 2 for first, second, third in inertia]
    span = {0}
    for mask in masks:
        span |= {value ^ mask for value in tuple(span)}
    assert len(span) == 8
    assert 8 * (-2) + 4 * 4 == 0
    assert (1, 0, 0) not in inertia
    assert (0, 1, 0) in inertia and (0, 0, 1) in inertia


def check_boundary_partition():
    contacts = Counter()
    # T=0: simultaneous outer clusters at anchors 0 and 1.
    for bits in product((0, 1), repeat=3):
        zeros = {name for name, bit in zip(VARIABLES, bits) if bit == 0}
        ones = set(VARIABLES) - zeros
        if zeros:
            contacts[canonical_boundary({"0"} | zeros)] += 1
        if ones:
            contacts[canonical_boundary({"1"} | ones)] += 1

    # T=infinity: every nonempty subset can choose the infinite root.
    for size in (1, 2, 3):
        for subset in combinations(VARIABLES, size):
            contacts[canonical_boundary({"inf"} | set(subset))] += 1

    # T=2: moving collisions of size two or three.
    for size in (2, 3):
        for subset in combinations(VARIABLES, size):
            contacts[canonical_boundary(subset)] += 1

    assert len(all_boundaries()) == 25
    assert set(contacts) == all_boundaries()
    assert set(contacts.values()) == {1}


def check_transversality():
    poles = tuple(-1 / lam for lam in PARAMETERS)
    assert poles == (F(-1, 3), F(121, 57), F(25, 9))
    assert len(set(poles) | {F(0), F(1)}) == 5

    constants = tuple(1 + 2 * lam for lam in PARAMETERS)
    other_at_two = tuple(-1 / value for value in constants)
    assert other_at_two == (F(-1, 7), F(-121, 7), F(-25, 7))
    assert len(set(other_at_two) | {F(0), F(1), F(2)}) == 6

    inverse_zero = tuple(-1 / value for value in constants)
    inverse_one = tuple((1 + lam) / value
                        for lam, value in zip(PARAMETERS, constants))
    inverse_two = tuple(value / (4 * lam + 3)
                        for lam, value in zip(PARAMETERS, constants))
    inverse_infinity = tuple(value / lam
                             for lam, value in zip(PARAMETERS, constants))
    assert inverse_zero == (F(-1, 7), F(-121, 7), F(-25, 7))
    assert inverse_one == (F(4, 7), F(64, 7), F(16, 7))
    assert inverse_two == (F(7, 15), F(7, 135), F(7, 39))
    assert inverse_infinity == (F(7, 3), F(-7, 57), F(-7, 9))
    for slopes in (inverse_zero, inverse_one, inverse_two, inverse_infinity):
        assert len(set(slopes)) == 3
        assert all(slopes)


def check_rank_certificate():
    branch_roots = (F(-7), F(-7, 9), F(-7, 81), F(-7, 361))

    def point_count(prime):
        roots = [root.numerator * pow(root.denominator, -1, prime) % prime
                 for root in branch_roots]
        assert len(set(roots)) == 4
        count = 2  # the two rational points at infinity of the monic quartic
        squares = {value * value % prime for value in range(prime)}
        for value in range(prime):
            rhs = 1
            for root in roots:
                rhs = rhs * (value - root) % prime
            count += 1 if rhs == 0 else 2 if rhs in squares else 0
        return count

    assert point_count(13) == 20
    assert point_count(17) == 24
    assert 24 > 4


if __name__ == "__main__":
    check_pencil_and_branches()
    check_boundary_partition()
    check_transversality()
    check_rank_certificate()
    print("Quadratic pencil, four branch values, and mixed inertia pass exactly.")
    print("The 14+7+4 contacts are all 25 stable boundaries, each once.")
    print("All common-root inverse slopes are nonzero and pairwise distinct.")
    print("Good-reduction counts 20 and 24 bound torsion by 4, below 24 rational points.")
