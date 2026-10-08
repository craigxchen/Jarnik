"""Exact covariance and scope checks for gauss_reduction_marked_nodes_audit.md."""

from fractions import Fraction as F
from functools import reduce
from itertools import product
from math import gcd, lcm, prod


def mv(m, v):
    return tuple(sum(a * b for a, b in zip(row, v)) for row in m)


def transpose(m):
    return tuple(zip(*m))


def mm(m, n):
    return tuple(tuple(sum(a * b for a, b in zip(row, col))
                       for col in transpose(n)) for row in m)


def det(v, w):
    return v[0] * w[1] - v[1] * w[0]


def inverse(s):
    (a, b), (c, d) = s
    assert a * d - b * c == 1
    return ((d, -b), (-c, a))


def qvalue(q, v):
    return sum(x * y for x, y in zip(v, mv(q, v)))


def weights(q, nodes):
    return [F(qvalue(q, v) ** 2,
              prod(det(v, w) for j, w in enumerate(nodes) if i != j))
            for i, v in enumerate(nodes)]


def primitive(values):
    scale = lcm(*(x.denominator for x in values))
    integers = [int(x * scale) for x in values]
    content = reduce(gcd, map(abs, integers))
    result = [x // content for x in integers]
    return result if result[0] > 0 else [-x for x in result]


def reduce_form(q):
    """Gauss shear/swap algorithm, keeping the exact cumulative SL2 matrix."""
    s = ((1, 0), (0, 1))
    for _ in range(10000):
        a, b, c = q[0][0], q[0][1], q[1][1]
        if abs(2 * b) > a:
            shift = -((2 * b + a) // (2 * a))
            step = ((1, shift), (0, 1))
        elif a > c:
            step = ((0, -1), (1, 0))
        else:
            return q, s
        q = mm(mm(transpose(step), q), step)
        s = mm(s, step)
    raise AssertionError("Reduction did not terminate")


def main():
    nodes = [(1, 0)] + [(t, 1) for t in (0, 1, 2, 3, 5)]
    q = ((2, 1), (1, 1))
    s = ((1, 0), (-1, 1))
    identity = ((1, 0), (0, 1))
    assert mm(mm(transpose(s), q), s) == identity
    moved = [mv(inverse(s), v) for v in nodes]
    assert [F(x, y) for x, y in moved] == list(map(F, (1, 0))) + [
        F(1, 2), F(2, 3), F(3, 4), F(5, 6)]
    expected = [-480, 4, -375, 3380, -6250, 3721]
    assert primitive(weights(q, nodes)) == [-x for x in expected]
    assert weights(q, nodes) == weights(identity, moved)

    count = 0
    for a, b, c, d in product(range(-4, 5), repeat=4):
        if a * d - b * c != 1:
            continue
        s = ((a, b), (c, d))
        q2 = mm(mm(transpose(s), q), s)
        moved = [mv(inverse(s), v) for v in nodes]
        assert weights(q, nodes) == weights(q2, moved)
        for v in moved:
            assert gcd(abs(v[0]), abs(v[1])) == 1
        reduced, cumulative = reduce_form(q2)
        assert reduced == identity
        assert mm(mm(transpose(cumulative), q2), cumulative) == reduced
        count += 1

    scales = list(map(F, (2, -3, 5, 7, 11, -13)))
    rescaled = [tuple(scale * x for x in v) for scale, v in zip(scales, nodes)]
    assert weights(q, rescaled) == [w / prod(scales) for w in weights(q, nodes)]
    assert primitive(weights(q, rescaled)) == primitive(weights(q, nodes))

    reduced_count = 0
    for m in range(1, 21):
        for a in range(1, 2 * m + 1):
            for b in range(-a // 2 - 1, a // 2 + 2):
                if (m * m + b * b) % a:
                    continue
                c = (m * m + b * b) // a
                if abs(2 * b) > a or c < a or gcd(gcd(a, b), c) != 1:
                    continue
                assert 3 * a * a <= 4 * m * m
                assert max(a, abs(b), c) <= 2 * m * m
                reduced_count += 1

    x = (1, 2, 3)
    assert all((x[i] * x[j] + 1) % (x[j] - x[i]) == 0
               for i in range(3) for j in range(i + 1, 3))
    assert F(x[1] - x[0], x[2] - x[0]) == F(1, 2)
    assert 5 * 5 - 4 * 4 == 9 and 9 % gcd(5, 5) != 0
    assert 4 * 5 - 4 * 4 == 4 and 2 % gcd(4, 4) != 0
    print(f"PASS: {count} exact SL2 covariance/reduction cases; {reduced_count} reduced forms")
    print("PASS: infinity-first vector, independent rescaling, and three small-determinant counterexamples")


if __name__ == "__main__":
    main()
