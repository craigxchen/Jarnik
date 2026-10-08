"""Exact checks for five_point_cm_norm_lift.md."""
from fractions import Fraction as F
from math import gcd
from random import Random
from itertools import combinations


def mul(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def conj(a):
    return a[0], -a[1]


def norm(a):
    return a[0] ** 2 + a[1] ** 2


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def lin(r, z, s, w):
    return r * z[0] + s * w[0], r * z[1] + s * w[1]


def bezout(a, b):
    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if old_r == -1:
        old_r, old_s, old_t = 1, -old_s, -old_t
    assert old_r == 1 and old_s * a + old_t * b == 1
    return old_s, old_t


def check_system(z, w, vs):
    assert det(z, w) == 1
    qs = [lin(r, z, s, w) for r, s in vs]
    assert all(det(qs[i], qs[j]) for i in range(5) for j in range(i))

    A, B, C = norm(z), dot(z, w), norm(w)
    assert mul(conj(z), w) == (B, 1)
    assert A * C - B * B == 1
    for (r, s), q in zip(vs, qs):
        X = dot(z, q)
        assert det(z, q) == s
        assert X == A * r + B * s
        assert (X - B * s) % A == 0
        assert X * X + s * s == A * norm(q)

    xs = [dot(z, q) for q in qs]
    ss = [v[1] for v in vs]
    for i in range(5):
        for j in range(i):
            assert xs[i] * ss[j] - ss[i] * xs[j] == A * det(qs[i], qs[j])

    d = lambda i, j: det(vs[i - 1], vs[j - 1])
    a = F(d(1, 3) * d(2, 4), d(1, 4) * d(2, 3))
    b = F(d(1, 3) * d(2, 5), d(1, 5) * d(2, 3))
    s2, s3, s4, s5 = (vs[j][1] for j in range(1, 5))
    for parameter, sj, qj in ((a, s4, qs[3]), (b, s5, qs[4])):
        left = tuple(F(s2 * s3 * x) for x in qj)
        right = tuple(
            sj * ((1 - parameter) * s3 * qs[1][k]
                  + parameter * s2 * qs[2][k])
            for k in range(2)
        )
        assert left == right

    assert a * s4 * (xs[1] * s3 - s2 * xs[2]) == \
        s3 * (xs[1] * s4 - s2 * xs[3])
    assert b * s5 * (xs[1] * s3 - s2 * xs[2]) == \
        s3 * (xs[1] * s5 - s2 * xs[4])


def check_contact_grouping():
    """Verify the median clusters and row divisors in (22)--(24)."""
    labels = set(range(1, 6))
    edges = list(combinations(range(1, 6), 2))
    clusters = {}
    row_edges = {i: set() for i in labels}
    for edge in edges:
        edge = set(edge)
        bits = [int(a in edge) for a in (1, 2, 3)]
        mu = sorted(bits)[1]
        cluster = {i for i in labels if int(i in edge) != mu}
        key = tuple(sorted(edge))
        clusters[key] = cluster
        for i in cluster:
            row_edges[i].add(key)

    expected_clusters = {
        (1, 4): {1, 4}, (1, 5): {1, 5}, (2, 3): {1, 4, 5},
        (2, 4): {2, 4}, (2, 5): {2, 5}, (1, 3): {2, 4, 5},
        (3, 4): {3, 4}, (3, 5): {3, 5}, (1, 2): {3, 4, 5},
        (4, 5): {4, 5},
    }
    assert clusters == expected_clusters
    assert row_edges == {
        1: {(1, 4), (1, 5), (2, 3)},
        2: {(2, 4), (2, 5), (1, 3)},
        3: {(3, 4), (3, 5), (1, 2)},
        4: {(1, 4), (2, 4), (3, 4), (2, 3), (1, 3), (1, 2), (4, 5)},
        5: {(1, 5), (2, 5), (3, 5), (2, 3), (1, 3), (1, 2), (4, 5)},
    }


def main():
    check_contact_grouping()
    rng = Random(20260912)
    checked = 0
    while checked < 2000:
        z = rng.randint(-30, 30), rng.randint(-30, 30)
        if gcd(*z) != 1:
            continue
        alpha, beta = bezout(*z)
        w0 = -beta, alpha
        shear = rng.randint(-20, 20)
        w = w0[0] + shear * z[0], w0[1] + shear * z[1]
        vs = [(1, 0)]
        while len(vs) < 5:
            v = rng.randint(-20, 20), rng.randint(-20, 20)
            if gcd(*v) != 1 or not v[1]:
                continue
            if any(det(v, old) == 0 for old in vs):
                continue
            vs.append(v)
        check_system(z, w, vs)
        checked += 1

    # Large congruence index alone does not force a long determinant-one pair.
    # For z=1,w=i and v_e=(Re G_e,Im G_e), v_e,1*z+v_e,2*w=G_e.
    gs = [(2, 1), (3, 2), (4, 1), (5, 2), (6, 1),
          (5, 4), (7, 2), (6, 5), (8, 3), (9, 4)]
    prime_norms = [5, 13, 17, 29, 37, 41, 53, 61, 73, 97]
    assert [norm(g) for g in gs] == prime_norms
    z, w = (1, 0), (0, 1)
    assert det(z, w) == 1
    assert all(lin(g[0], z, g[1], w) == g for g in gs)
    assert all(gcd(*g) == 1 for g in gs)

    print(f"{checked} exact determinant-one norm-conic and five-point lifts pass.")
    print("The ten median contact clusters group as three, three, three, seven, seven row factors.")
    print("A ten-congruence lattice can have arbitrarily large index and a unit-size determinant-one solution.")


if __name__ == "__main__":
    main()
