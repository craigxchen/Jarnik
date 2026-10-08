"""Finite checks for six_row_elliptic_forgetful_obstruction.md."""
from fractions import Fraction
from itertools import combinations, product


def canonical_boundary(part, labels):
    part = frozenset(part)
    other = frozenset(labels) - part
    return min(tuple(sorted(part)), tuple(sorted(other)))


def check_boundaries():
    labels = set(range(6))
    boundaries = {
        canonical_boundary(part, labels)
        for size in (2, 3, 4)
        for part in combinations(labels, size)
        if 2 <= len(labels - set(part))
    }
    assert len(boundaries) == 25

    for forgotten in labels:
        kept = labels - {forgotten}
        downstairs = list(combinations(sorted(kept), 2))
        assert len(downstairs) == 10
        for edge in downstairs:
            lifts = {
                canonical_boundary(edge, labels),
                canonical_boundary(set(edge) | {forgotten}, labels),
            }
            assert len(lifts) == 2
            assert lifts <= boundaries

    # A boundary point of M_0,4 pulls back under forgetting two labels to
    # the four ways of distributing those labels across its 2+2 split.
    for retained_tuple in combinations(sorted(labels), 4):
        retained = set(retained_tuple)
        forgotten = sorted(labels - retained)
        first_pair = set(retained_tuple[:2])
        fiber = {
            canonical_boundary(first_pair | {
                forgotten[index]
                for index in range(2)
                if mask & (1 << index)
            }, labels)
            for mask in range(4)
        }
        assert len(fiber) == 4
        assert fiber <= boundaries


def check_del_pezzo_class():
    # Intersections with E_i give m_i=delta.  The six joining lines then
    # give degree-2delta=delta, hence degree=3delta.
    for delta in range(1, 20):
        multiplicities = [delta] * 4
        degree = 3 * delta
        assert all(degree - multiplicities[i] - multiplicities[j] == delta
                   for i, j in combinations(range(4), 2))
        genus = 1 + 5 * delta * (delta - 1) // 2
        if delta == 1:
            assert genus == 1
        else:
            assert genus > 1

    # Exact d=1 alternatives and their singularity budgets.
    assert 1 + 5 * 1 * 0 // 2 - 1 == 0  # smooth anticanonical
    assert 1 + 5 * 1 * 0 // 2 - 0 == 1  # rational anticanonical
    assert 1 + 5 * 2 * 1 // 2 - 1 == 5  # genus-one |2L| normalization


def check_kernel_intersection():
    # Nonzero E[2] points are bit vectors 1,2,3.  Model only the selected
    # two-torsion contained in each even-order kernel; enlarging kernels
    # can only enlarge the pair sums and their intersection.
    nonzero = (1, 2, 3)
    for ta, tb, tc in product(nonzero, repeat=3):
        ka, kb, kc = {0, ta}, {0, tb}, {0, tc}
        add = lambda first, second: {x ^ y for x in first for y in second}
        intersection = add(kb, kc) & add(ka, kc) & add(ka, kb)
        assert intersection - {0}


def check_involution_caps():
    # Nonzero vectors of F_2^3 model involutions in a rank-three
    # elementary abelian subgroup.  Three distinct nonzero vectors are
    # dependent exactly when their xor is zero.
    vectors = range(1, 8)
    caps = []
    for size in range(1, 8):
        for subset in combinations(vectors, size):
            if all(first ^ second ^ third
                   for first, second, third in combinations(subset, 3)):
                caps.append(subset)
    assert max(map(len, caps)) == 4
    four_caps = [subset for subset in caps if len(subset) == 4]
    assert all(first ^ second ^ third ^ fourth == 0
               for first, second, third, fourth in four_caps)

    # Take {0,1,2,3} as the translation plane and {4,5,6,7} as its
    # reflection coset.  If every pair quotient is rational, no pair in
    # the cap can consist of translations.  The unique size-four option
    # is then the full reflection coset.
    translations = {1, 2, 3}
    rational_pair_caps = [
        subset for subset in four_caps
        if sum(value in translations for value in subset) <= 1
    ]
    assert rational_pair_caps == [(4, 5, 6, 7)]


def check_reflection_cap_boundary_collision():
    def sigma(value):
        return value / (2 * value - 1)

    # The nontrivial common projectivity is an involution, has parameter
    # -1, and fixes precisely 0 and 1.
    for value in map(Fraction, (-3, -1, 0, 1, 2, 4)):
        if 2 * value == 1:
            continue
        image = sigma(value)
        if 2 * image != 1:
            assert sigma(image) == value
        if value not in (0, 1):
            assert image != value
        if value not in (0, 1) and image not in (0, 1):
            parameter = image * (1 - value) / (value * (1 - image))
            assert parameter == -1

    labels = {"inf", "0", "1", "a", "b", "c"}
    for pattern in product((0, 1), repeat=3):
        zero_coordinates = {
            label for label, value in zip(("a", "b", "c"), pattern)
            if value == 0
        }
        if zero_coordinates:
            boundary_side = {"0"} | zero_coordinates
        else:
            boundary_side = {"1", "a", "b", "c"}
        assert 2 <= len(boundary_side) <= 4
        assert 2 <= len(labels - boundary_side) <= 4


def check_three_reflection_quadratics():
    def f_value(A, B, C, alpha, x):
        v = (A - C) * alpha * alpha
        u = (B - C) * (1 - alpha) ** 2 - v
        return C + (u * x + v) / (x - alpha) ** 2

    samples = tuple(map(Fraction, (-3, -2, -1, 2, 3, 4)))
    for A, B, C in product(map(Fraction, (-1, 0, 2)), repeat=3):
        if A == B == C:
            continue
        for alpha in samples:
            assert f_value(A, B, C, alpha, Fraction(0)) == A
            assert f_value(A, B, C, alpha, Fraction(1)) == B

            # The derivative at zero is a nonconstant quadratic in
            # 1/alpha unless all three anchor values agree.
            derivative_zero = (
                (B - A) + 2 * (A - B) / alpha
                + (B - C) / (alpha * alpha)
            )
            v = (A - C) * alpha * alpha
            u = (B - C) * (1 - alpha) ** 2 - v
            derivative_direct = u / (alpha * alpha) + 2 * v / alpha**3
            assert derivative_direct == derivative_zero

        assert (2 * (A - B), B - C) != (0, 0)
        assert (-2 * (B - C), B - A) != (0, 0)

        # For a nonanchor x, the quadratic f_alpha(x)=y cannot vanish
        # identically in alpha unless A=B=C.  Exhaust the small exact grid.
        for x in map(Fraction, (-2, -1, 2, 3)):
            for y in map(Fraction, (-2, -1, 0, 1, 2)):
                coefficients = (
                    (C - y) * x * x + x * (B - C),
                    -2 * x * (B - y),
                    (C - y) + x * (B - C) + (A - C) * (1 - x),
                )
                assert coefficients != (0, 0, 0)


def main():
    check_boundaries()
    check_del_pezzo_class()
    check_kernel_intersection()
    check_involution_caps()
    check_reflection_cap_boundary_collision()
    check_three_reflection_quadratics()
    print("All 25 boundaries and every two-term forgetful pullback pass.")
    print("Equal ten-line intersections force the stated anticanonical multiple.")
    print("Every triple of even elliptic kernels has nontrivial pair-sum intersection.")
    print("Degree-one deck involutions form caps of size at most four.")
    print("A four-reflection cap forces four copies of one boundary collision.")
    print("Three common-branch reflection quotients fail the quadratic collision test.")


if __name__ == "__main__":
    main()
