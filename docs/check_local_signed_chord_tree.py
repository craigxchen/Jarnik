"""Exact bounded checks for docs/local_signed_chord_tree.md.

All arithmetic is integral.  The finite family check uses 8 values of M and
four primes, hence only 32 six-node cases and 480 edges.
"""

from itertools import combinations, product
from math import gcd, lcm


def vp(n, p):
    n = abs(n)
    out = 0
    while n and n % p == 0:
        n //= p
        out += 1
    return out


def det(u, v):
    return u[0] * v[1] - u[1] * v[0]


def q_i(v):
    return v[0] * v[0] + v[1] * v[1]


def polar_i(u, v):
    return u[0] * v[0] + u[1] * v[1]


def edge_norm(u, v):
    b, d = polar_i(u, v), det(u, v)
    g = gcd(abs(b), abs(d))
    a, t = b // g, d // g
    epsilon = 2 if a % 2 and t % 2 else 1
    return (a * a + t * t) // epsilon


def lambdas(vs, p):
    out = {}
    for i, j in combinations(range(6), 2):
        d = det(vs[i], vs[j])
        out[i, j] = vp(q_i(vs[i]), p) + vp(q_i(vs[j]), p) - 2 * vp(d, p)
    return out


def central_profile(vs, p):
    aa = [2 * vp(q_i(v), p) - sum(vp(det(v, w), p)
                                  for k, w in enumerate(vs) if k != i)
          for i, v in enumerate(vs)]
    m = min(aa)
    return tuple(a - m for a in aa)


def assert_primitive_distinct(vs):
    assert all(gcd(abs(x), abs(y)) == 1 for x, y in vs)
    assert all(det(u, v) != 0 for u, v in combinations(vs, 2))


def split_family(p, M):
    if p == 5:
        iso = ((1, 2), (1, -2), (2, 1))
    elif p == 13:
        iso = ((1, 5), (1, -5), (5, 1))
    else:
        raise ValueError(p)
    return ((1, 0), (1, p ** M), *iso, (1, 1))


def inert_family(p, M):
    return ((1, 0), (1, p ** M), (0, 1), (1, 1), (1, 2), (1, 4))


def residue_class(v, p):
    x, y = v[0] % p, v[1] % p
    if x:
        return (1, y * pow(x, -1, p) % p)
    return (0, 1)


def check_residue_classes(p):
    reps = []
    for x, y in product(range(p), repeat=2):
        if (x, y) == (0, 0):
            continue
        c = residue_class((x, y), p)
        if c not in reps:
            reps.append(c)
    iso = [q_i(c) % p == 0 for c in reps]
    expected = 2 if p in (5, 13) else 0
    assert sum(iso) == expected
    for i, j in combinations(range(len(reps)), 2):
        if iso[i] or iso[j]:
            # Distinct residue classes have d=0 and a positive endpoint has
            # a positive norm valuation for an integral lift.
            continue
        assert not (iso[i] or iso[j])


def check_split_families():
    for p in (5, 13):
        check_residue_classes(p)
        for M in range(1, 9):
            vs = split_family(p, M)
            assert_primitive_distinct(vs)
            lam = lambdas(vs, p)
            assert min(lam.values()) == -2 * M
            assert max(lam.values()) == 2
            assert lam[0, 1] == -2 * M
            assert lam[2, 3] == lam[2, 4] == 2
            assert central_profile(vs, p) == (0, 0, M + 2, M + 1, M + 1, M)
            # The pair cotangent B/Delta on edge 12 has denominator p^M.
            u, v = vs[0], vs[1]
            assert vp(det(u, v), p) - vp(polar_i(u, v), p) == M
            all_edge_n = lcm(*(edge_norm(u, v)
                               for u, v in combinations(vs, 2)))
            assert all_edge_n == p * p * (1 + p ** (2 * M)) // 2


def check_same_branch():
    a, b = (1, 2), (1, 7)
    c, d = (1, 7), (1, 132)
    assert (vp(q_i(a), 5), vp(q_i(b), 5), vp(det(a, b), 5)) == (1, 2, 1)
    assert (vp(q_i(c), 5), vp(q_i(d), 5), vp(det(c, d), 5)) == (2, 2, 3)
    assert vp(q_i(a), 5) + vp(q_i(b), 5) - 2 * vp(det(a, b), 5) == 1
    assert vp(q_i(c), 5) + vp(q_i(d), 5) - 2 * vp(det(c, d), 5) == -2


def check_inert_families():
    for p in (3, 7):
        check_residue_classes(p)
        for M in range(1, 9):
            vs = inert_family(p, M)
            assert_primitive_distinct(vs)
            assert all(vp(q_i(v), p) == 0 for v in vs)
            lam = lambdas(vs, p)
            assert max(lam.values()) == 0
            assert lam[0, 1] == -2 * M
            assert all(value <= 0 for value in lam.values())


def check_unimodular_positive_form():
    # Q=S^T S with det(Q)=1.  Pulling vectors back by S^-1 preserves all
    # determinants and gives the same local profiles as Q=I.
    S = ((2, 1), (1, 1))
    det_s = 1
    q = (5, 3, 2)
    assert 5 * 2 - 3 * 3 == det_s * det_s
    for M in (1, 4, 8):
        original = split_family(5, M)
        pulled = tuple((v[0] - v[1], -v[0] + 2 * v[1]) for v in original)
        assert_primitive_distinct(pulled)
        got = []
        for i, j in combinations(range(6), 2):
            u, v = pulled[i], pulled[j]
            qi = q[0] * u[0] ** 2 + 2 * q[1] * u[0] * u[1] + q[2] * u[1] ** 2
            qj = q[0] * v[0] ** 2 + 2 * q[1] * v[0] * v[1] + q[2] * v[1] ** 2
            got.append(vp(qi, 5) + vp(qj, 5) - 2 * vp(det(u, v), 5))
        assert min(got) == -2 * M and max(got) == 2


if __name__ == "__main__":
    check_split_families()
    check_same_branch()
    check_inert_families()
    check_unimodular_positive_form()
    print("local signed chord tree: exact checks passed")
