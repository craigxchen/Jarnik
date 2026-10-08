"""Exact whole-distance Smith invariants on actual Gaussian circle tuples."""

from functools import reduce
from itertools import combinations
from math import gcd, isqrt, prod
from random import Random

from check_quartet_matching_gcd_cut_budget import ggcd, mul, norm, sub


def conj(z):
    return z[0], -z[1]


def gp(z, n):
    out = (1, 0)
    for _ in range(n):
        out = mul(out, z)
    return out


def divide(z, w):
    numerator, denominator = mul(z, conj(w)), norm(w)
    assert all(x % denominator == 0 for x in numerator)
    return tuple(x // denominator for x in numerator)


def primitive(points):
    common = reduce(ggcd, points)
    return [divide(z, common) for z in points]


def det2(a):
    return a[0][0] * a[1][1] - a[0][1] * a[1][0]


def det3(a):
    return sum((-1) ** j * a[0][j]
               * det2([[a[r][c] for c in range(3) if c != j]
                       for r in (1, 2)]) for j in range(3))


def minor(E, I, J):
    matrix = [[E[i][j] for j in J] for i in I]
    return det2(matrix) if len(I) == 2 else det3(matrix)


def audit(points, common_unit=False):
    m, N = len(points), norm(points[0])
    assert m >= 3 and len(set(points)) == m
    assert all(norm(z) == N for z in points)
    E = [[N - z[0] * w[0] - z[1] * w[1] for w in points]
         for z in points]
    pairs, triples = list(combinations(range(m), 2)), list(combinations(range(m), 3))
    D = {}
    for i, j, k in triples:
        u, v = sub(points[j], points[i]), sub(points[k], points[i])
        D[i, j, k] = u[0] * v[1] - u[1] * v[0]
        assert D[i, j, k]
    d = gcd(*(E[i][j] for i, j in pairs))
    g = gcd(*D.values())
    delta2 = gcd(*(minor(E, I, J) for I in pairs for J in pairs))
    delta3 = gcd(*(minor(E, I, J) for I in triples for J in triples))
    assert delta2 == d * d
    assert delta3 == N * g * g
    assert delta3 % (d ** 3) == 0
    assert g % d == 0
    for I in triples:
        for J in triples:
            assert minor(E, I, J) == N * D[I] * D[J]
    is_primitive = norm(reduce(ggcd, points)) == 1
    if is_primitive:
        assert gcd(N, d) == 1 and (g * g) % (d ** 3) == 0
        parity_types = {(x % 2, y % 2) for x, y in points}
        kappa = 1 if len(parity_types) == 2 else 2
        assert d % kappa == 0
        h_free = isqrt(d // kappa)
        assert d == kappa * h_free * h_free
        assert g % (2 * kappa * h_free ** 3) == 0
        a_free = g // (2 * kappa * h_free ** 3)
        assert delta3 // d ** 3 == (4 // kappa) * N * a_free * a_free
    if not common_unit:
        return d, g, None
    assert is_primitive and N % 2
    residues = {}
    for i, j in pairs:
        n = norm(ggcd(points[i], points[j]))
        assert E[i][j] % (2 * n) == 0
        square = E[i][j] // (2 * n)
        t = isqrt(square)
        assert t > 0 and t * t == square
        residues[i, j] = t
    h = gcd(*residues.values())
    assert d == 2 * h * h and g % (4 * h ** 3) == 0
    a = g // (4 * h ** 3)
    for I in triples:
        n = norm(reduce(ggcd, (points[i] for i in I)))
        assert abs(D[I]) == 4 * n * prod(residues[e] for e in combinations(I, 2))
    assert delta3 // d ** 3 == 2 * N * a * a
    return d, g, h


def main():
    rng = Random(20261001)
    split = [(2, 1), (3, 2), (4, 1), (5, 2)]
    units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    common_count = mixed_count = scaled_count = 0
    for case in range(24):
        width = 2 + case % 3
        exponents = [1 + rng.randrange(4) for _ in range(width)]
        m = 3 + case % 6
        possibilities = prod(e + 1 for e in exponents)
        m = min(m, possibilities)
        allocations = {tuple(0 for _ in exponents), tuple(exponents)}
        while len(allocations) < m:
            allocations.add(tuple(rng.randrange(e + 1) for e in exponents))
        rows = []
        for allocation in sorted(allocations):
            z = (1, 0)
            for pi, e, b in zip(split, exponents, allocation):
                z = mul(z, mul(gp(pi, b), gp(conj(pi), e - b)))
            rows.append(z)
        epsilon = units[case % 4]
        rows = [mul(epsilon, z) for z in rows]
        audit(rows, True)
        common_count += 1
        mixed = [mul(units[rng.randrange(4)], z) for z in rows]
        assert len(set(mixed)) == len(mixed)
        audit(mixed)
        mixed_count += 1
        audit([mul((2 + case % 2, 1), z) for z in mixed])
        scaled_count += 1

    half_angles = [(x, 6) for x in (1, 5, 7, 23, 35, 37, 41)]
    product = reduce(mul, half_angles, (1, 0))
    rows = [conj(product)] + [mul(divide(conj(product), conj(z)), z) for z in half_angles]
    rows = primitive(rows)
    _, _, h = audit(rows, True)
    assert h == 6
    mixed_half_angles = [(1, 3), (2, 3), (4, 3)]
    product = reduce(mul, mixed_half_angles, (1, 0))
    rows = primitive([conj(product)]
                     + [mul(divide(conj(product), conj(z)), z) for z in mixed_half_angles])
    assert audit(rows)[0] == 9
    # A primitive mixed-parity circle need not have the common-unit factor two.
    assert audit([(1, 0), (0, 1), (-1, 0), (0, -1)])[0] == 1
    print(f"verified {common_count} nested common-unit, {mixed_count} mixed-unit, "
          f"{scaled_count} scaled tuples; h=6 and mixed-parity h=3 fixtures")


if __name__ == "__main__":
    main()
