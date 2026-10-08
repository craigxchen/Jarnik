"""Exact audit of joint Gaussian gcd normalization on the Pell fixture."""

from itertools import combinations, product
from math import gcd, isqrt, lcm, prod

from check_gaussian_near_real_divisor_reduction import pair_factor
from check_gaussian_reflection_replacement import (
    gconj, gdivexact, gcd_many, gmul, gnorm, ggcd,
)


BLOCKS = [(-1, 2), (13, 8), (8, 5), (50, 31)]
ROWS = [(1, 1, 1, 1), (0, 1, 0, 0),
        (0, 0, 1, 0), (0, 0, 0, 1)]
EDGES = list(combinations(range(len(ROWS)), 2))


def content(z):
    return gcd(abs(z[0]), abs(z[1]))


def lift(a, n):
    d = gnorm(a)
    assert n % d == 0
    aa = gmul(a, a)
    return (n // d * aa[0], n // d * aa[1])


def fixture_data():
    n = prod(gnorm(block) for block in BLOCKS)
    assert n == 358853785
    factors = {}
    lifts = {}
    for edge in EDGES:
        a = pair_factor(BLOCKS, ROWS[edge[0]], ROWS[edge[1]])
        assert gcd(a[0], a[1]) == 1
        assert gnorm(ggcd(a, gconj(a))) == 1
        factors[edge] = a
        lifts[edge] = lift(a, n)
        assert gnorm(lifts[edge]) == n * n
    return n, factors, lifts


def check_ordinary_content_identity(n, lifts, factors):
    checked = 0
    source = []
    for row in ROWS:
        z = (1, 0)
        for block, exponent in zip(BLOCKS, row):
            z = gmul(z, block if exponent else gconj(block))
        source.append(z)
    assert all(gnorm(z) == n for z in source)
    for size in (2, 3, 4):
        for subset in combinations(range(len(ROWS)), size):
            family = [edge for edge in EDGES
                      if edge[0] in subset and edge[1] in subset]
            integer_content = gcd(*(content(lifts[edge]) for edge in family))
            pair_lcm = lcm(*(gnorm(factors[edge]) for edge in family))
            assert n % pair_lcm == 0
            assert integer_content == n // pair_lcm
            common = gcd_many(source[i] for i in subset)
            assert integer_content == gnorm(common)
            reduced = {i: gdivexact(source[i], common) for i in subset}
            for i,j in family:
                assert lifts[i,j] == gmul(source[i], gconj(source[j]))
                assert gdivexact(lifts[i,j], (integer_content,0)) == gmul(reduced[i], gconj(reduced[j]))
            checked += 1
    return checked


def check_positive_three_point_case(n, lifts):
    family = [(0, 1), (0, 2), (1, 2)]
    positive = [z if z[1] > 0 else gconj(z) for z in
                (lifts[edge] for edge in family)]
    assert gcd(*(content(z) for z in positive)) == 1
    assert gcd(*(v for z in positive for v in z)) == 1
    gaussian = gcd_many(positive)
    assert gaussian == (-144, -1)
    assert gnorm(gaussian) == 20737
    normalized = [gdivexact(z, gaussian) for z in positive]
    normalized_norm = n * n // gnorm(gaussian)
    assert normalized_norm == 6209964749425
    assert all(gnorm(z) == normalized_norm for z in normalized)
    assert isqrt(normalized_norm) ** 2 != normalized_norm
    assert n * n % gnorm(gaussian) == 0
    return gaussian, normalized_norm


def check_integer_radius_oriented_case(n, factors, lifts):
    # The three edges on source rows (1,2,3) have a mixed orientation:
    # Z_12, conjugate(Z_13), Z_23.  It has a rational gcd but does not
    # preserve a two-factor product relation without conjugation.
    edges = [(1, 2), (1, 3), (2, 3)]
    oriented = [lifts[edges[0]], gconj(lifts[edges[1]]), lifts[edges[2]]]
    gaussian = gcd_many(oriented)
    assert gaussian == (-5, 0)
    divisor = (5, 0)  # an associate of the computed gcd
    assert all(gdivexact(z, divisor) for z in oriented)
    radius = n // 5
    normalized = [gdivexact(z, divisor) for z in oriented]
    assert radius == 71770757
    assert all(gnorm(z) == radius * radius for z in normalized)
    for i, j, k in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
        product = gmul(normalized[i], normalized[j])
        assert product not in (
            (radius * normalized[k][0], radius * normalized[k][1]),
            (-radius * normalized[k][0], -radius * normalized[k][1]),
        )
    for z in normalized:
        h = radius - z[0]
        assert h > 0
        assert z[1] * z[1] == h * (2 * radius - h)
    assert [radius - z[0] for z in normalized] == [6922, 1602, 1864]
    # Cyclic orientation still retains an exact three-factor relation.
    assert gmul(gmul(normalized[0],normalized[1]),normalized[2]) == (radius**3,0)
    return gaussian, radius


def relation_with_coefficient(values, coefficient):
    for i, j, k in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
        product = gmul(values[i], values[j])
        target = (coefficient * values[k][0], coefficient * values[k][1])
        if product == target or product == (-target[0], -target[1]):
            return True
    return False


def check_all_orientations(n, lifts):
    edges = [(1, 2), (1, 3), (2, 3)]
    raw = [lifts[edge] for edge in edges]
    square_norm_with_relation = []
    for signs in product((1, -1), repeat=3):
        values = [z if sign == 1 else gconj(z)
                  for z, sign in zip(raw, signs)]
        common = gcd_many(values)
        if relation_with_coefficient(values, n):
            square_norm_with_relation.append(gnorm(common))
    assert square_norm_with_relation
    assert all(isqrt(norm) ** 2 != norm for norm in square_norm_with_relation)
    return len(square_norm_with_relation)


def main():
    n, factors, lifts = fixture_data()
    identities = check_ordinary_content_identity(n, lifts, factors)
    nonsquare_gcd, nonsquare_norm = check_positive_three_point_case(n, lifts)
    integer_gcd, integer_radius = check_integer_radius_oriented_case(
        n, factors, lifts)
    relation_orientations = check_all_orientations(n, lifts)
    print(f"PASS: N={n}, {len(EDGES)} Pell pair lifts, all norm N^2.")
    print(f"PASS: ordinary content identity on {identities} source subsets.")
    print(f"PASS: positive-imaginary subset gcd {nonsquare_gcd}, norm "
          f"{gnorm(nonsquare_gcd)}, normalized norm {nonsquare_norm}.")
    print(f"PASS: oriented rational gcd {integer_gcd}, integer radius "
          f"{integer_radius}, and three exact individual Pell equations.")
    print(f"PASS: all 8 orientations checked; {relation_orientations} retain "
          "a product relation, all with nonsquare gcd norms.")


if __name__ == "__main__":
    main()
