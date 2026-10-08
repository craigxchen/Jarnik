"""Exact integer audit of the moving-quadrilateral fifth-point identities."""

from itertools import combinations
from math import atan2, gcd, prod
from random import Random


def integer_determinant(rows):
    matrix = [list(row) for row in rows]
    previous, sign = 1, 1
    for pivot in range(len(matrix) - 1):
        other = next(i for i in range(pivot, len(matrix)) if matrix[i][pivot])
        if other != pivot:
            matrix[pivot], matrix[other] = matrix[other], matrix[pivot]
            sign = -sign
        value = matrix[pivot][pivot]
        for i in range(pivot + 1, len(matrix)):
            for j in range(pivot + 1, len(matrix)):
                numerator = value * matrix[i][j] - matrix[i][pivot] * matrix[pivot][j]
                assert numerator % previous == 0
                matrix[i][j] = numerator // previous
            matrix[i][pivot] = 0
        previous = value
    return sign * matrix[-1][-1]


def audit(points, check_determinant=False):
    def det(i, j, k):
        u = tuple(points[j][h] - points[i][h] for h in range(2))
        v = tuple(points[k][h] - points[i][h] for h in range(2))
        return u[0] * v[1] - u[1] * v[0]

    def chord2(i, j):
        return sum((points[i][h] - points[j][h]) ** 2 for h in range(2))

    norm = sum(x * x for x in points[0])
    A, B, C, E = det(0, 1, 2), det(0, 1, 3), det(0, 2, 3), det(1, 2, 3)
    assert min(A, B, C, E) > 0 and E == A - B + C
    pi = A * B * C * E
    M = gcd(B * chord2(0, 2), C * chord2(0, 1))
    assert A * B * C % M == 0
    a, b, s = B * chord2(0, 2) // M, C * chord2(0, 1) // M, A * B * C // M
    assert a > b > 0 and s > 0 and gcd(a, b) == 1
    K = B * C + A * E
    assert K * K - (B - A) ** 2 * (A + C) ** 2 == 4 * pi
    assert -(B - A) ** 2 * a * a + 2 * K * a * b - (A + C) ** 2 * b * b == 4 * s * s
    assert 4 * s ** 3 * norm == pi * a * b * (a - b)
    assert s * chord2(1, 3) == a * B * E
    F, G = det(0, 1, 4) * det(2, 3, 4), det(0, 2, 4) * det(1, 3, 4)
    k = gcd(abs(F), abs(G))
    assert k > 0 and a * F == b * G
    assert abs(F) == b * k and abs(G) == a * k
    assert pi * prod(chord2(i, 4) for i in range(4)) == 16 * k * k * s * s * norm * norm
    if check_determinant:
        rows = [(1, x, y, x * x, x * y) for x, y in points]
        assert abs(integer_determinant(rows)) == k * s


def main():
    circles = {}
    for x in range(-44, 45):
        for y in range(-44, 45):
            norm = x * x + y * y
            if 0 < norm <= 2000:
                circles.setdefault(norm, []).append((x, y))
    rng = Random(20260914)
    count = 0
    for norm, points in sorted(circles.items()):
        if len(points) < 5:
            continue
        points.sort(key=lambda p: atan2(p[1], p[0]))
        samples = combinations(points, 5) if len(points) <= 12 else (
            sorted(rng.sample(points, 5), key=lambda p: atan2(p[1], p[0]))
            for _ in range(300)
        )
        for sample in samples:
            # Each of the five choices for the extra point checks its sign
            # relative to four cyclically ordered base points.
            for extra in range(5):
                audit([p for i, p in enumerate(sample) if i != extra] + [sample[extra]],
                      check_determinant=(count % 997 == 0))
                count += 1
    print(f"Verified {count:,} ordered quadrilateral/fifth-point configurations.")


if __name__ == "__main__":
    main()
