"""Exact Farkas certificates for five of the seven balanced M=6 cuts.

This is a finite combinatorial check.  It does not classify other positive
dependencies, odd numbers of rows, or multilevel allocations.
"""

from itertools import combinations, permutations, product


M = 6
DEPENDENCE = (1, 1, 1, 1, 2, 2)
PAIRS = tuple(combinations(range(M), 2))


def balanced_classes():
    """Return quotient representatives of cuts of dependence mass 4."""
    out = []
    for mask in range(1, 1 << M):
        cut = tuple(i for i in range(M) if mask & (1 << i))
        if sum(DEPENDENCE[i] for i in cut) != 4:
            continue
        complement = tuple(i for i in range(M) if not mask & (1 << i))
        if cut <= complement:
            out.append(cut)
    return tuple(out)


def sigma(pair, cut):
    i, j = pair
    return 1 if ((i in cut) != (j in cut)) else -1


def certificate(code):
    """Find y>=0 with y^T sigma <= 0 and a strict negative coefficient."""
    matrix = tuple(
        tuple(sigma(pair, cut) for cut in code)
        for pair in PAIRS
    )
    # First search small positive integer coefficients for a certificate
    # negative in every column.  The found coefficients are at most 3.  This
    # distinguishes weak infeasibility from a system whose weak solutions lie
    # on a face.  Then search 0/1 coefficients for the latter certificates.
    for support_size, require_all_negative, value_vectors in (
        (5, True, set(permutations((1, 1, 2, 2, 3)))),
        (1, False, ((1,),)), (2, False, ((1, 1),)),
        (3, False, ((1, 1, 1),)),
    ):
        for support in combinations(range(len(PAIRS)), support_size):
            for values in value_vectors:
                coefficients = tuple(
                    sum(value * matrix[q][column]
                        for q, value in zip(support, values))
                    for column in range(len(code))
                )
                if max(coefficients) <= 0 and min(coefficients) < 0 and (
                    not require_all_negative or max(coefficients) < 0
                ):
                    return support, values, coefficients
    # This line is unreachable for the present seven classes, but keeps the
    # function's failure mode explicit if the finite instance is edited.
    raise AssertionError("no certificate found")


def main():
    classes = balanced_classes()
    assert classes == (
        (0, 1, 2, 3),
        (0, 1, 4), (0, 2, 4), (0, 3, 4),
        (0, 1, 5), (0, 2, 5), (0, 3, 5),
    )

    # The complete profile has a strict witness: all seven weights are 1.
    full_margins = tuple(
        sum(sigma(pair, cut) for cut in classes) for pair in PAIRS
    )
    assert set(full_margins) == {1, 5}
    assert all(value > 0 for value in full_margins)

    counts = {"weakly_infeasible": 0, "boundary_only": 0}
    for code in combinations(classes, 5):
        support, values, coefficients = certificate(code)
        assert all(value > 0 for value in values)
        assert tuple(sum(value * sigma(PAIRS[q], cut)
                         for q, value in zip(support, values))
                     for cut in code) == coefficients
        assert all(value <= 0 for value in coefficients)
        assert any(value < 0 for value in coefficients)
        if all(value < 0 for value in coefficients):
            counts["weakly_infeasible"] += 1
        else:
            assert coefficients.count(0) == 4
            assert coefficients.count(-2) == 1
            counts["boundary_only"] += 1
            # A certificate alone does not establish existence on the boundary.
            # Verify a nonzero weak witness separately; divide by four to
            # normalize its total weight to one.
            assert any(
                sum(weights) == 4 and
                all(sum(sigma(pair, cut) * weight
                        for cut, weight in zip(code, weights)) >= 0
                    for pair in PAIRS)
                for weights in product((0, 1), repeat=5)
            )

        # Exact contradiction: strict margins and a nonnegative certificate
        # make the same sum positive, while its displayed coefficient form is
        # strictly negative for every strictly positive column-weight vector.
        assert sum(coefficients) < 0

    assert counts == {"weakly_infeasible": 12, "boundary_only": 9}
    print(
        "PASS: 7 balanced classes; all 21 five-column subsets have exact "
        "Farkas certificates (12 weakly infeasible, 9 boundary-only); "
        "full seven-column weights (1,...,1) are strictly obtuse."
    )


if __name__ == "__main__":
    main()
