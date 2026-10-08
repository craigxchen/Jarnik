"""Exact finite checks for four_row_nearest_square_gram_rigidity.md."""

from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt, prod


def dot(x, y):
    return x[0] * y[0] + x[1] * y[1]


def det(x, y):
    return x[0] * y[1] - x[1] * y[0]


def norm(x):
    return dot(x, x)


def mul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def gauss_quotient(x, y):
    n = norm(y)
    z = mul(x, (y[0], -y[1]))
    return (z[0] // n, z[1] // n) if z[0] % n == z[1] % n == 0 else None


def gauss_gcd(x, y):
    while y != (0, 0):
        n = norm(y)
        z = mul(x, (y[0], -y[1]))
        q = (round(Fraction(z[0], n)), round(Fraction(z[1], n)))
        x, y = y, (x[0] - mul(q, y)[0], x[1] - mul(q, y)[1])
    return x


def gcd_many(points):
    ans = (0, 0)
    for point in points:
        ans = gauss_gcd(ans, point)
    return ans


def matrix_rank(rows):
    a = [[Fraction(x) for x in row] for row in rows]
    rank = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(rank, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][col]
        a[rank] = [x / scale for x in a[rank]]
        for i in range(len(a)):
            if i != rank:
                scale = a[i][col]
                a[i] = [x - scale * y for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def determinant(a):
    if len(a) == 1:
        return a[0][0]
    return sum((-1) ** j * a[0][j] * determinant(
        [row[:j] + row[j + 1:] for row in a[1:]])
        for j in range(len(a)))


def cofactor_relation(b):
    """Return a nonzero small integer vector killed by rank-2/3 B."""
    rank = matrix_rank(b)
    assert rank in (2, 3)
    rows = next(rs for rs in combinations(range(4), rank)
                if matrix_rank([b[i] for i in rs]) == rank)
    for cols in combinations(range(4), rank + 1):
        c = [0] * 4
        for k, col in enumerate(cols):
            sub = [[b[i][j] for j in cols if j != col] for i in rows]
            c[col] = (-1) ** k * determinant(sub)
        if any(c):
            assert all(sum(b[i][j] * c[j] for j in range(4)) == 0
                       for i in range(4))
            return c
    raise AssertionError("missing cofactor relation")


def gram_data(points):
    n = [norm(p) for p in points]
    a = [isqrt(v) for v in n]
    b = [[dot(points[i], points[j]) - a[i] * a[j]
          for j in range(4)] for i in range(4)]
    return n, a, b


def relation_threshold(points, b):
    g = gcd_many(points)
    reduced = [gauss_quotient(p, g) for p in points]
    assert all(p is not None for p in reduced)
    d = [gcd_many([reduced[j] for j in range(4) if j != i])
         for i in range(4)]
    e = max(1, *(abs(x) for row in b for x in row))
    l_squared = min(norm(x) for x in d)
    return (6 * e ** 3) ** 2 < l_squared, e, l_squared


def check_cofactor_barrier():
    # Both higher ranks really occur for positive integral Gram matrices.
    fixtures = [
        ([(2, 0), (2, 3), (3, 0), (3, 3)], 2),
        ([(4, 1), (4, 3), (5, 4), (8, 1)], 3),
    ]
    for points, expected_rank in fixtures:
        _, _, b = gram_data(points)
        assert matrix_rank(b) == expected_rank
        c = cofactor_relation(b)
        assert all(sum(c[i] * points[i][j] for i in range(4)) == 0
                   for j in (0, 1))
        triggers, e, l_squared = relation_threshold(points, b)
        assert not triggers
        assert max(map(abs, c)) <= 6 * e ** 3
        assert max(x * x for x in c) >= l_squared


def check_integral_rotation():
    # Q is near the real axis; P=(3-4i)Q/5 is integral and has the same
    # oriented Gram data.  This tests the rational rotation in the proof.
    q = [(1002, 1), (1004, 2), (1011, 3), (1023, 4)]
    p = []
    for point in q:
        product = mul((3, -4), point)
        assert product[0] % 5 == product[1] % 5 == 0
        p.append((product[0] // 5, product[1] // 5))
    n, a, b = gram_data(p)
    assert a == [x for x, _ in q]
    assert [n[i] - a[i] ** 2 for i in range(4)] == [y * y for _, y in q]
    assert matrix_rank(b) == 1
    assert all(b[i][j] == q[i][1] * q[j][1]
               for i in range(4) for j in range(4))
    assert all(dot(p[i], p[j]) == dot(q[i], q[j])
               and det(p[i], p[j]) == det(q[i], q[j])
               for i, j in combinations(range(4), 2))
    # Recover the rational unit vector u from a=H^t u.
    d = det(p[0], p[1])
    assert d
    u1 = Fraction(a[0] * p[1][1] - a[1] * p[0][1], d)
    u2 = Fraction(p[0][0] * a[1] - p[1][0] * a[0], d)
    assert u1 * u1 + u2 * u2 == 1
    assert all(u1 * x + u2 * y == a[i]
               and -u2 * x + u1 * y == q[i][1]
               for i, (x, y) in enumerate(p))
    for j in range(1, 4):
        z = (dot(p[0], p[j]), det(p[0], p[j]))
        assert gauss_quotient(z, (a[0], -q[0][1])) == q[j]


def crt(congruences):
    modulus = prod(m for _, m in congruences)
    value = sum((r % m) * (modulus // m) * pow(modulus // m, -1, m)
                for r, m in congruences)
    return value % modulus, modulus


def check_triggered_threshold():
    # Four co-singleton Gaussian primes make the finite inequality true.
    # The common rational unit rotation keeps the Gram data integral.
    factors = [(5, 4), (7, 2), (5, 6), (8, 3)]
    prime_norms = [norm(f) for f in factors]
    assert prime_norms == [41, 53, 61, 73]
    roots = [f[0] * pow(f[1], -1, norm(f)) % norm(f) for f in factors]
    q = []
    for j in range(4):
        conditions = [(2, 5), (0, 2)]
        conditions += [((root + (i == j)) % prime, prime)
                       for i, (root, prime) in enumerate(zip(roots, prime_norms))]
        x, modulus = crt(conditions)
        q.append((x + (j + 1) * modulus, 1))
    assert all(gauss_quotient(q[j], factors[i]) is not None
               for i in range(4) for j in range(4) if i != j)
    assert all(gauss_quotient(q[i], factors[i]) is None for i in range(4))
    p = []
    for point in q:
        product = mul((3, -4), point)
        assert product[0] % 5 == product[1] % 5 == 0
        p.append((product[0] // 5, product[1] // 5))
    _, a, b = gram_data(p)
    assert a == [x for x, _ in q]
    assert all(x == 1 for row in b for x in row)
    triggers, e, l_squared = relation_threshold(p, b)
    assert triggers and e == 1 and l_squared > 36
    assert all(det(p[i], p[j]) == det(q[i], q[j])
               for i, j in combinations(range(4), 2))


def check_prime_power_trimming():
    # One full-core prime survives after trimming eight of ten powers;
    # all fifteen distinct block supports remain represented.
    primes = [5, 13, 17, 29, 37, 41, 53, 61, 73, 89,
              97, 101, 109, 113, 137]
    subsets = [frozenset(s) for size in range(1, 5)
               for s in combinations(range(4), size)]
    core = {t: p ** (10 if len(t) == 4 else 1)
            for t, p in zip(subsets, primes)}
    # The last assigned prime belongs to the full block.
    full_prime = primes[-1]
    content = [full_prime] * 4
    c = prod(content)
    trimmed = {t: value // gcd(value, c * c) for t, value in core.items()}
    assert trimmed[frozenset(range(4))] == full_prime ** 2
    for i in range(4):
        old_norm = prod(value for t, value in core.items() if i in t)
        new_norm = old_norm // content[i] ** 2
        assert old_norm % content[i] ** 2 == 0
        assert all(new_norm % value == 0 for t, value in trimmed.items() if i in t)
    for i, j in combinations(range(4), 2):
        old_delta = prod(value for t, value in core.items() if i in t and j in t)
        new_delta = old_delta // (content[i] * content[j])
        assert old_delta % (content[i] * content[j]) == 0
        assert all(new_delta % value == 0 for t, value in trimmed.items()
                   if i in t and j in t)


def check_three_column_pell_limit():
    seeds = [(4, 1, 2), (4, 3, 2), (5, 4, 2)]
    expected = [[13, 15, 20], [15, 21, 28], [20, 28, 37]]
    assert determinant(expected) == -16
    assert matrix_rank(expected) == 3
    u, v = 3, 2
    for _ in range(12):
        assert u * u - 2 * v * v == 1
        j = [[(u + 1) // 2, (u - 1) // 2, v],
             [(u - 1) // 2, (u + 1) // 2, v], [v, v, u]]
        lifted = [tuple(sum(row[k] * seed[k] for k in range(3))
                        for row in j) for seed in seeds]
        points = [p[:2] for p in lifted]
        a = [p[2] for p in lifted]
        assert a == [5 * v + 2 * u, 7 * v + 2 * u, 9 * v + 2 * u]
        assert [isqrt(norm(p)) for p in points] == a
        assert [[dot(points[i], points[k]) - a[i] * a[k]
                 for k in range(3)] for i in range(3)] == expected
        assert [norm(points[i]) - a[i] ** 2 for i in range(3)] == [13, 21, 37]
        minors = [det(points[i], points[k]) for i, k in combinations(range(3), 2)]
        assert minors == [8 * u + 4 * v, 11 * u + 4 * v, u]
        assert gcd(gcd(*minors[:2]), minors[2]) == 1
        relation = [u, -11 * u - 4 * v, 8 * u + 4 * v]
        assert all(sum(relation[i] * points[i][k] for i in range(3)) == 0
                   for k in range(2))
        u, v = 3 * u + 4 * v, 2 * u + 3 * v


if __name__ == "__main__":
    check_cofactor_barrier()
    check_integral_rotation()
    check_triggered_threshold()
    check_prime_power_trimming()
    check_three_column_pell_limit()
    print("PASS: cofactor barrier, triggered threshold, integral rotation, prime-power trimming, and Pell limit")
