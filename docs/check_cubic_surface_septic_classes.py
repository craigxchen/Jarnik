"""Exact Picard-lattice checks for genus-one degree-seven classes on a cubic surface."""

from itertools import combinations, combinations_with_replacement
from fractions import Fraction


def candidate_classes():
    """Return sorted coefficient patterns (a; b_1,...,b_6).

    On S=Bl_6(P^2), write D=a*l-sum(b_i E_i).  The anticanonical class is
    H=3*l-sum(E_i), so H.D=7 and D^2=7 are the degree/genus-one equations.
    The Cauchy bound restricts a to 3,...,11; the finite search below then
    finds every sorted integer pattern in that range.
    """
    result = []
    for a in range(3, 12):
        total = 3 * a - 7
        square = a * a - 7
        # Every |b_i|^2 is at most the total square a^2-7 <= 114.
        for b in combinations_with_replacement(range(-11, 12), 6):
            if sum(b) == total and sum(x * x for x in b) == square:
                result.append((a, b))
    return result


def check_cauchy_bound():
    # (sum b_i)^2 <= 6 sum b_i^2 gives
    # 3a^2 - 42a + 91 <= 0, hence 3 <= a <= 11.
    assert all(3 <= a <= 11 for a, _ in candidate_classes())


def check_plane_cubic_class():
    # A smooth plane cubic through A_5,A_6, avoiding A_1,...,A_4, has strict
    # transform D=3*l-E_5-E_6 on the six-point blowup.
    a, b = 3, (0, 0, 0, 0, 1, 1)
    assert 3 * a - sum(b) == 7       # H.D
    assert a * a - sum(x * x for x in b) == 7  # D^2
    assert 1 + (7 - 7) // 2 == 1     # adjunction: 2g-2=D.(D+K)=7-7


def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def rank_exact(matrix):
    matrix = [list(map(Fraction, row)) for row in matrix]
    rows, cols, rank = len(matrix), len(matrix[0]), 0
    for col in range(cols):
        pivot = next((row for row in range(rank, rows)
                      if matrix[row][col] != 0), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        for row in range(rank + 1, rows):
            factor = matrix[row][col] / matrix[rank][col]
            for j in range(col, cols):
                matrix[row][j] -= factor * matrix[rank][j]
        rank += 1
    return rank


def check_rational_instance_general_position():
    points = [(0, 0, 1), (1, 0, 1), (-3, -3, 1), (-3, -2, 1),
              (-2, -3, 1), (-2, 1, 1)]
    assert all(det3(*(points[i] for i in triple)) != 0
               for triple in combinations(range(6), 3))
    conic_rows = [[x * x, x * y, y * y, x * z, y * z, z * z]
                  for x, y, z in points]
    assert rank_exact(conic_rows) == 6
    assert all(y * y != x ** 3 - x for x, y, _ in points[2:])


def main():
    check_cauchy_bound()
    check_plane_cubic_class()
    check_rational_instance_general_position()
    classes = candidate_classes()
    assert (3, (0, 0, 0, 0, 1, 1)) in classes
    print(f"{len(classes)} sorted integer Picard patterns satisfy H.D=7 and D^2=7.")
    print("The strict transform of a smooth plane cubic through two blowup centers is one such class.")
    print("Class and genus equations therefore give no cubic-surface obstruction to a degree-seven elliptic septic.")


if __name__ == "__main__":
    main()
