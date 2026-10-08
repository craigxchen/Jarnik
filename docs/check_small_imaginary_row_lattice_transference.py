"""Exact checks supplementing small_imaginary_row_lattice_transference.md."""

from itertools import combinations
from math import gcd, prod


def norm(z):
    return z[0] ** 2 + z[1] ** 2


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def det(z, w):
    return z[0] * w[1] - z[1] * w[0]


def quotient(z, w):
    n = norm(w)
    x, y = mul(z, (w[0], -w[1]))
    assert x % n == y % n == 0
    return x // n, y // n


def ggcd(z, w):
    while w != (0, 0):
        n = norm(w)
        x, y = mul(z, (w[0], -w[1]))
        q = ((2 * x + n) // (2 * n), (2 * y + n) // (2 * n))
        qw = mul(q, w)
        z, w = w, (z[0] - qw[0], z[1] - qw[1])
    return z


def primitive(z):
    return gcd(*z) == 1 and (z[0] + z[1]) % 2 == 1


def check_tuple(rows, blocks):
    assert all(primitive(z) for z in rows)
    common = (0, 0)
    for z in rows:
        common = ggcd(common, z)
    d = norm(common)
    qs = [quotient(z, common) for z in rows]
    pairs = list(combinations(range(len(rows)), 2))
    ds = [det(rows[i], rows[j]) for i, j in pairs]
    assert all(ds)
    f = gcd(*ds)
    assert f % d == 0
    h = f // d
    assert h == gcd(*(det(qs[i], qs[j]) for i, j in pairs))
    assert gcd(h, prod(map(norm, qs))) == 1
    ts = []
    for (i, j), delta in zip(pairs, ds):
        core = prod(n for subset, n in blocks.items() if i in subset and j in subset)
        assert delta % core == 0
        ts.append(delta // core)
    assert all(t % h == 0 for t in ts)
    gy = gcd(*(z[1] for z in rows))
    assert h % gy == 0
    return d, h, ts


def check_family():
    for k in range(1, 41):
        a = 2210 * k
        h0 = (a, 1)
        qs = [(2 * a + 1, -2), (4 * a + 1, -4)]
        d = norm(h0)
        ns = [norm(q) for q in qs]
        assert all(gcd(x, y) == 1 for x, y in combinations([d] + ns, 2))
        assert all(primitive(z) for z in [h0] + qs)
        rows = [mul(h0, q) for q in qs]
        assert rows == [(2 * d + a, 1), (4 * d + a, 1)]
        blocks = {frozenset([0, 1]): d, frozenset([0]): ns[0], frozenset([1]): ns[1]}
        assert check_tuple(rows, blocks) == (d, 2, [-2])
        assert (2 * qs[0][0] - qs[1][0], 2 * qs[0][1] - qs[1][1]) == (1, 0)
        c1, c2 = 4 * a + 1, -(2 * a + 1)
        assert (c1 * qs[0][0] + c2 * qs[1][0], c1 * qs[0][1] + c2 * qs[1][1]) == (0, 2)
        assert abs(c1) + abs(c2) == 6 * a + 2
    return 40


def check_three_row_fixture():
    rows = [(1469178, 1), (31686, 1), (153548, 1)]
    blocks = {
        frozenset([0, 1, 2]): 109,
        frozenset([0, 1]): 157,
        frozenset([1, 2]): 13,
        frozenset([0, 2]): 85,
        frozenset([0]): 1483897,
        frozenset([1]): 4513,
        frozenset([2]): 195749,
    }
    assert all(gcd(x, y) == 1 for x, y in combinations(blocks.values(), 2))
    for i, z in enumerate(rows):
        assert norm(z) == prod(n for subset, n in blocks.items() if i in subset)
    assert check_tuple(rows, blocks) == (109, 2, [84, 142, -86])


def check_nonisotropic_index():
    rows = [(x, y) for x in range(-4, 5) for y in range(-4, 5) if primitive((x, y))]
    checked = 0
    for size in (2, 3):
        for zs in combinations(rows, size):
            ds = [det(z, w) for z, w in combinations(zs, 2)]
            if not all(ds):
                continue
            common = (0, 0)
            for z in zs:
                common = ggcd(common, z)
            if norm(common) != 1:
                continue
            h = gcd(*ds)
            assert gcd(h, prod(map(norm, zs))) == 1
            checked += 1
    return checked


if __name__ == "__main__":
    count = check_family()
    check_three_row_fixture()
    small = check_nonisotropic_index()
    print(f"Passed {count} balanced family instances, one full-support three-row fixture, and {small} primitive index checks.")
