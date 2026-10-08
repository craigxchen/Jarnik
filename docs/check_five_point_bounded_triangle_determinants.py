"""Exact checks for five-point determinant rigidity and its near-circle obstruction."""

from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import atan2, gcd, isqrt


def determinant(matrix):
    size = len(matrix)
    if size == 0:
        return 1
    return sum(
        (-1) ** column * matrix[0][column]
        * determinant([row[:column] + row[column + 1:] for row in matrix[1:]])
        for column in range(size)
    )


def det(left, right):
    return left[0] * right[1] - left[1] * right[0]


def sub(left, right):
    return (left[0] - right[0], left[1] - right[1])


def dot(left, right):
    return left[0] * right[0] + left[1] * right[1]


def check_five_points():
    norm = 1105
    nodes = []
    for x in range(1, isqrt(norm) + 1):
        y = isqrt(norm - x * x)
        if y > 0 and x * x + y * y == norm:
            nodes.append((x, y))
    nodes.sort(key=lambda point: atan2(point[1], point[0]))
    checked = 0
    for points in combinations(nodes, 5):
        edges = [sub(points[i + 1], points[i]) for i in range(4)]
        adjacent = [det(edges[i], edges[i + 1]) for i in range(3)]
        assert all(value > 0 for value in adjacent)
        ceiling = max(adjacent)
        bound = 2 ** 18 * ceiling ** 10
        assert all(abs(det(a, b)) <= bound for a, b in combinations(edges, 2))
        d = adjacent[0]
        marks = []
        for point in points:
            difference = sub(point, points[0])
            marks.append((det(difference, edges[1]), det(edges[0], difference)))
        assert all(abs(coordinate) <= 4 * bound for mark in marks for coordinate in mark)
        matrix = [[x * x, x * y, y * y, x, y] for x, y in marks[1:]]
        kernel = [
            (-1) ** column * determinant([row[:column] + row[column + 1:] for row in matrix])
            for column in range(5)
        ]
        content = reduce(gcd, kernel)
        assert content != 0
        kernel = [value // abs(content) for value in kernel]
        if kernel[0] < 0:
            kernel = [-value for value in kernel]
        a, b, c, u, v = map(Fraction, kernel)
        coefficient_bound = 24 * (16 * bound ** 2) ** 4
        assert all(abs(value) <= coefficient_bound for value in kernel)
        discriminant = a * c - b * b / 4
        assert a > 0 and discriminant >= Fraction(1, 4)
        scale = Fraction(dot(edges[0], edges[0]), d * d) / a
        assert scale * b / 2 == Fraction(dot(edges[0], edges[1]), d * d)
        assert scale * c == Fraction(dot(edges[1], edges[1]), d * d)
        assert scale * scale * discriminant == Fraction(1, d * d)
        assert scale <= Fraction(2, d)
        radius_squared = scale * (c * u * u - b * u * v + a * v * v) / (4 * discriminant)
        assert radius_squared == norm
        assert radius_squared <= Fraction(8 * coefficient_bound ** 3, d)
        checked += 1
    return checked


def check_near_circle():
    checked = 0
    for parameter in (3, 11, 101, 10001):
        count = isqrt(isqrt(parameter))
        center = (
            (parameter * parameter + 1) // 2,
            -parameter * (parameter * parameter + 1) // 2,
        )
        radius_squared = (parameter * parameter + 1) ** 3 // 4
        points = [
            (parameter * j + j * j - center[0], j - center[1])
            for j in range(count + 1)
        ]
        for j, point in enumerate(points):
            defect = j ** 3 * (2 * parameter + j)
            assert dot(point, point) - radius_squared == defect
            assert defect <= 3 * parameter * count ** 3
            assert dot(points[0], point) > 0
        edges = [sub(points[j + 1], points[j]) for j in range(count)]
        for j, edge in enumerate(edges):
            assert edge == (parameter + 2 * j + 1, 1)
            assert gcd(*edge) == 1
        for j in range(count - 1):
            assert det(edges[j], edges[j + 1]) == -2
        for j in range(1, count - 1):
            assert edges[j + 1] == tuple(2 * edges[j][axis] - edges[j - 1][axis] for axis in range(2))
        checked += len(points)
    return checked


if __name__ == "__main__":
    tuples = check_five_points()
    near_points = check_near_circle()
    print("PASS: {} exact five-point metric reconstructions, including determinant 1/d^2 and radius scale.".format(tuples))
    print("PASS: {} integer near-circle points with exact norm defects, primitive chords, and bounded recurrences.".format(near_points))
