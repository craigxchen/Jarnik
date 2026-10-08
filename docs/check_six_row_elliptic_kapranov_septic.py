"""Combinatorial checks for six_row_elliptic_kapranov_septic.md."""
from itertools import combinations


FOUR = frozenset(range(4))


def p(index):
    return ("p", index)


def q(first, second):
    return ("q",) + tuple(sorted((first, second)))


def r(*indices):
    return ("r",) + tuple(sorted(indices))


def coordinate_divisor(index):
    complement = sorted(FOUR - {index})
    return (
        {p(other) for other in complement}
        | {q(first, second) for first, second in combinations(complement, 2)}
        | {r(*complement)}
    )


def check_coordinate_gcds():
    divisors = [coordinate_divisor(index) for index in FOUR]
    assert all(len(divisor) == 7 for divisor in divisors)
    for first, second in combinations(FOUR, 2):
        other_first, other_second = sorted(FOUR - {first, second})
        forced_gcd = (
            divisors[other_first] & divisors[other_second]
        ) - {r(*sorted(FOUR - {other_first})),
             r(*sorted(FOUR - {other_second}))}
        assert forced_gcd == {p(first), p(second), q(first, second)}


def check_ten_plane_sections():
    planes = set()
    # x_i=0 contains the other three coordinate centers.
    for index in FOUR:
        planes.add(tuple(sorted(FOUR - {index})))
    # x_i-x_j=0 contains P_5 and the two complementary centers.
    for first, second in combinations(FOUR, 2):
        planes.add(tuple(sorted((FOUR - {first, second}) | {4})))
    assert len(planes) == 10
    assert planes == set(combinations(range(5), 3))

    # The six forced zeros of x_i-x_j are exactly the three centers and
    # the three line contacts in that plane; the seventh zero is residual.
    for first, second in combinations(FOUR, 2):
        third, fourth = sorted(FOUR - {first, second})
        forced = {
            p(4), p(third), p(fourth),
            q(4, third), q(4, fourth), q(third, fourth),
        }
        divisor = forced | {r(4, third, fourth)}
        assert len(divisor) == 7


def check_quadric_frame_lines():
    # A symmetric 4-by-4 coefficient array has ten independent entries.
    # Containment of every coordinate-frame line kills the two diagonal
    # and one off-diagonal coefficient indexed by each pair.
    killed = set()
    for first, second in combinations(FOUR, 2):
        killed |= {(first, first), (second, second), (first, second)}
    coefficients = {(i, j) for i in FOUR for j in FOUR if i <= j}
    assert killed == coefficients


if __name__ == "__main__":
    check_coordinate_gcds()
    check_ten_plane_sections()
    check_quadric_frame_lines()
    print("Four coordinate divisors force the six frame-line trisecants.")
    print("Four coordinate and six difference sections give all ten planes.")
    print("The six lines among four frame points kill every quadric coefficient.")
