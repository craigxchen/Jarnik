"""Exact fixtures for the four-witness matching-norm kernel bound.

The checker exercises the general determinant identities on retained exact
Gaussian-coordinate samples, then checks the matching factors and quantitative
inequality on signed Gaussian-product fixtures.
"""

from itertools import combinations, permutations
from math import gcd
import random


Gaussian = tuple[int, int]


def gmul(z: Gaussian, w: Gaussian) -> Gaussian:
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def gconj(z: Gaussian) -> Gaussian:
    return z[0], -z[1]


def gnorm(z: Gaussian) -> int:
    return z[0] * z[0] + z[1] * z[1]


def determinant(u: Gaussian, v: Gaussian) -> int:
    return u[0] * v[1] - v[0] * u[1]


def deltas(points: tuple[Gaussian, ...]) -> dict[tuple[int, int], int]:
    return {(i, j): determinant(points[i], points[j])
            for i, j in combinations(range(4), 2)}


def triple_rows(points: tuple[Gaussian, ...], t: dict[str, int]) -> tuple[tuple[int, int, int], ...]:
    y = tuple(z[1] for z in points)
    return (
        (y[2] * t["01"], -y[1] * t["02"], y[0] * t["12"]),
        (y[3] * t["01"], y[0] * t["13"], -y[1] * t["03"]),
        (y[0] * t["23"], y[3] * t["02"], -y[2] * t["03"]),
        (y[1] * t["23"], -y[2] * t["13"], y[3] * t["12"]),
    )


def dot3(u: tuple[int, int, int], v: tuple[int, int, int]) -> int:
    return sum(a * b for a, b in zip(u, v))


def cross3(u: tuple[int, int, int], v: tuple[int, int, int]) -> tuple[int, int, int]:
    return (u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0])


def matching_t(d: dict[tuple[int, int], int], scale: int,
               a: int, b: int, c: int) -> dict[str, int]:
    denominators = {"01": scale * a, "23": scale * a,
                    "02": scale * b, "13": scale * b,
                    "03": scale * c, "12": scale * c}
    out = {}
    for key, den in denominators.items():
        value = d[tuple(map(int, key))]
        assert den > 0 and value % den == 0
        out[key] = value // den
        assert out[key] != 0
    return out


def check_matching(points: tuple[Gaussian, ...], d0: int,
                   a: int, b: int, c: int) -> tuple[int, int]:
    """Check six paired divisibilities, triple rows, rank, and height bound."""
    assert all(z[1] != 0 for z in points)
    delta = deltas(points)
    assert all(value != 0 for value in delta.values())
    assert d0 > 0 and a > 0 and b > 0 and c > 0
    assert gcd(a, b) == gcd(a, c) == gcd(b, c) == 1
    t = matching_t(delta, d0, a, b, c)
    pairs = (("01", "23", a), ("02", "13", b), ("03", "12", c))
    for first, second, factor in pairs:
        assert delta[tuple(map(int, first))] == d0 * factor * t[first]
        assert delta[tuple(map(int, second))] == d0 * factor * t[second]

    rows = triple_rows(points, t)
    kernel = (a, b, c)
    assert all(dot3(row, kernel) == 0 for row in rows)
    cross_products = [cross3(rows[i], rows[j])
                      for i, j in combinations(range(4), 2)]
    nonzero_crosses = [v for v in cross_products if v != (0, 0, 0)]
    assert nonzero_crosses  # rank is at least two
    assert all(v[0] * b == v[1] * a and v[0] * c == v[2] * a
               for v in nonzero_crosses)
    # Pairwise coprimality makes (a,b,c) primitive; every integer cross
    # product parallel to it must therefore be an integer multiple of it.
    for v in nonzero_crosses:
        assert v[0] % a == 0
        k = v[0] // a
        assert v == (k * a, k * b, k * c)
        assert k != 0

    y_max = max(abs(z[1]) for z in points)
    t_max = max(abs(v) for v in t.values())
    assert max(a, b, c) <= 2 * (y_max * t_max) ** 2
    return len(nonzero_crosses), t_max


def factor_integer(n: int) -> dict[int, int]:
    n = abs(n)
    factors = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            factors[p] = factors.get(p, 0) + 1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def pairwise_coprime_factors(
        pair_gcds: tuple[int, int, int]) -> tuple[int, tuple[int, int, int]]:
    """Choose pairwise-coprime positive factors, one from each pair gcd."""
    common = gcd(gcd(pair_gcds[0], pair_gcds[1]), pair_gcds[2])
    residuals = tuple(g // common for g in pair_gcds)
    factors = [factor_integer(n) for n in residuals]
    chosen = [1, 1, 1]
    primes = sorted(set().union(*(f.keys() for f in factors)))
    for p in primes:
        exponents = [f.get(p, 0) for f in factors]
        owner = max(range(3), key=lambda i: (exponents[i], -i))
        if exponents[owner]:
            chosen[owner] *= p ** exponents[owner]
    assert gcd(chosen[0], chosen[1]) == gcd(chosen[0], chosen[2]) == gcd(chosen[1], chosen[2]) == 1
    assert all(residuals[i] % chosen[i] == 0 for i in range(3))
    return common, tuple(chosen)


def check_random_coordinate_fixtures() -> tuple[int, int, int]:
    rng = random.Random(20261006)
    retained = nontrivial = cross_count = 0
    for _ in range(12000):
        points = tuple((rng.randint(-12, 12), rng.choice(tuple(range(-12, 0)) + tuple(range(1, 13))))
                       for _ in range(4))
        delta = deltas(points)
        if any(v == 0 for v in delta.values()):
            continue
        pair_gcds = (gcd(abs(delta[(0, 1)]), abs(delta[(2, 3)])),
                     gcd(abs(delta[(0, 2)]), abs(delta[(1, 3)])),
                     gcd(abs(delta[(0, 3)]), abs(delta[(1, 2)])))
        scale, (a, b, c) = pairwise_coprime_factors(pair_gcds)
        rank_pairs, _ = check_matching(points, scale, a, b, c)
        retained += 1
        nontrivial += (a > 1 or b > 1 or c > 1)
        cross_count += rank_pairs
        if retained == 500:
            break
    assert retained == 500
    return retained, nontrivial, cross_count


def signed_products(k: int, d: Gaussian, a: Gaussian, b: Gaussian,
                    c: Gaussian) -> tuple[Gaussian, ...]:
    factors = ((a, b, c), (a, gconj(b), gconj(c)),
               (gconj(a), b, gconj(c)), (gconj(a), gconj(b), c))
    out = []
    for row in factors:
        z = (k * d[0], k * d[1])
        for factor in row:
            z = gmul(z, factor)
        out.append(z)
    return tuple(out)


def signed_products_with_corrections(corrections: tuple[Gaussian, ...], d: Gaussian,
                                    a: Gaussian, b: Gaussian,
                                    c: Gaussian) -> tuple[Gaussian, ...]:
    factors = ((a, b, c), (a, gconj(b), gconj(c)),
               (gconj(a), b, gconj(c)), (gconj(a), gconj(b), c))
    out = []
    for correction, row in zip(corrections, factors):
        z = gmul(correction, d)
        for factor in row:
            z = gmul(z, factor)
        out.append(z)
    return tuple(out)


def check_varying_correction_fixture(corrections: tuple[Gaussian, ...],
                                     d: Gaussian, a: Gaussian,
                                     b: Gaussian, c: Gaussian) -> tuple[bool, bool, bool]:
    """Check the strong and universal bounds for one Gaussian-product quartet."""
    points = signed_products_with_corrections(corrections, d, a, b, c)
    n_a, n_b, n_c = map(gnorm, (a, b, c))
    n_d = gnorm(d)
    m_squared = max(map(gnorm, corrections))
    y_max = max(abs(z[1]) for z in points)
    universal_rhs = 8 * m_squared * max(1, y_max) ** 4
    assert min(n_a, n_b, n_c) <= universal_rhs

    delta = deltas(points)
    strong_case = all(z[1] != 0 for z in points) and all(v != 0 for v in delta.values())
    max_norm_branch = strong_case and min(n_a, n_b, n_c) > m_squared
    if strong_case:
        check_matching(points, n_d, n_a, n_b, n_c)
        assert n_d <= 8 * m_squared * y_max ** 4
        if max_norm_branch:
            assert max(map(gnorm, points)) <= 8 * m_squared * y_max ** 4
    return strong_case, any(v == 0 for v in delta.values()), max_norm_branch


def check_varying_corrections_and_exceptions() -> tuple[int, int, int, int, int]:
    candidates = ((2, 1), (3, 2), (4, 1), (5, 2),
                  (4, 3), (6, 1), (5, 4), (7, 2))
    correction_pool = ((1, 0), (0, 1), (1, 1), (1, -1),
                       (2, 0), (0, 2), (2, 1), (2, -1))
    rng = random.Random(8122026)
    tested = strong = d_unit = determinant_zero = max_norm_branch = 0
    for a, b, c in combinations(candidates, 3):
        na, nb, nc = map(gnorm, (a, b, c))
        if gcd(na, nb) != 1 or gcd(na, nc) != 1 or gcd(nb, nc) != 1:
            continue
        d_choices = [(1, 0)] + [z for z in candidates
                                if all(gcd(gnorm(z), n) == 1 for n in (na, nb, nc))]
        for d in d_choices:
            for _ in range(16):
                corrections = tuple(rng.choice(correction_pool) for _ in range(4))
                if len(set(corrections)) == 1:
                    continue
                has_strong, has_zero_delta, has_max_norm = check_varying_correction_fixture(
                    corrections, d, a, b, c)
                tested += 1
                strong += has_strong
                d_unit += (gnorm(d) == 1 and has_strong)
                determinant_zero += has_zero_delta
                max_norm_branch += has_max_norm

    # Force a zero imaginary coordinate: K_0=conj(D*S_0) makes W_0 real.
    a, b, c, d = (2, 1), (3, 2), (4, 1), (1, 0)
    base = signed_products(1, d, a, b, c)
    real_correction = gconj(base[0])
    corrections = (real_correction, (1, 0), (1, 1), (1, -1))
    points = signed_products_with_corrections(corrections, d, a, b, c)
    assert points[0][1] == 0
    check_varying_correction_fixture(corrections, d, a, b, c)

    # Force Δ_01=0 by choosing K_0=S_1 and K_1=S_0, so W_0=W_1.
    signed = ((a, b, c), (a, gconj(b), gconj(c)),
              (gconj(a), b, gconj(c)), (gconj(a), gconj(b), c))
    def product3(row):
        z = (1, 0)
        for factor in row:
            z = gmul(z, factor)
        return z

    corrections = (product3(signed[1]), product3(signed[0]), (1, 1), (1, -1))
    points = signed_products_with_corrections(corrections, d, a, b, c)
    assert points[0] == points[1]
    assert deltas(points)[(0, 1)] == 0
    check_varying_correction_fixture(corrections, d, a, b, c)
    return tested, strong, d_unit, determinant_zero, max_norm_branch, 2


def check_exact_group_reorientations() -> int:
    """Exhaust all column orders/signs and choices of the common group."""
    base = ((1, 1, 1, 1),
            (1, 1, -1, -1),
            (1, -1, 1, -1),
            (1, -1, -1, 1))
    checks = 0
    for order in permutations(range(4)):
        for column_mask in range(16):
            source = tuple(tuple(row[order[j]] *
                                 (-1 if (column_mask >> j) & 1 else 1)
                                 for j in range(4)) for row in base)
            for common_column in range(4):
                oriented = source
                column = tuple(row[common_column] for row in oriented)
                # A whole-column conjugation chooses the desired orientation
                # for this group. Otherwise it is already balanced and the
                # two negative rows are conjugated to make it common.
                if all(x == -1 for x in column):
                    oriented = tuple(tuple(-v if j == common_column else v
                                            for j, v in enumerate(row))
                                     for row in oriented)
                    column = tuple(row[common_column] for row in oriented)
                flipped_rows = tuple(i for i, row in enumerate(oriented)
                                     if row[common_column] == -1)
                assert len(flipped_rows) in (0, 2)
                normalized = tuple(tuple(-v for v in row) if i in flipped_rows else row
                                   for i, row in enumerate(oriented))
                assert all(row[common_column] == 1 for row in normalized)
                assert all(sum(x * y for x, y in zip(u, v)) == 0
                           for u, v in combinations(normalized, 2))
                varying = tuple(j for j in range(4) if j != common_column)
                patterns = {tuple(row[j] for j in varying) for row in normalized}
                assert len(patterns) == 4
                row_parities = {p[0] * p[1] * p[2] for p in patterns}
                assert len(row_parities) == 1
                for col in varying:
                    assert sum(row[col] == 1 for row in normalized) == 2
                checks += 1
    return checks


def check_gaussian_product_fixtures() -> tuple[int, int]:
    candidates = ((2, 1), (3, 2), (4, 1), (5, 2),
                  (4, 3), (6, 1), (5, 4), (7, 2))
    accepted = nonzero_bound = 0
    for a_block, b_block, c_block in combinations(candidates, 3):
        na, nb, nc = map(gnorm, (a_block, b_block, c_block))
        if gcd(na, nb) != 1 or gcd(na, nc) != 1 or gcd(nb, nc) != 1:
            continue
        for d_block in candidates:
            nd = gnorm(d_block)
            if any(gcd(nd, n) != 1 for n in (na, nb, nc)):
                continue
            for k in (1, 2, 3):
                points = signed_products(k, d_block, a_block, b_block, c_block)
                if any(z[1] == 0 for z in points):
                    continue
                delta = deltas(points)
                if any(v == 0 for v in delta.values()):
                    continue
                scale = k * k * nd
                t = matching_t(delta, scale, na, nb, nc)
                # These are the exact group-norm divisibilities, not factors
                # selected after seeing the pair gcds.
                check_matching(points, scale, na, nb, nc)
                y_max = max(abs(z[1]) for z in points)
                n_plus, n_minus = max(na, nb, nc), min(na, nb, nc)
                # Here all four corrections are the same positive integer k,
                # so M=k; nD is the norm of the displayed common block D.
                if k < n_minus:
                    lhs_scaled = k * k * y_max ** 4 * (8 * na * nb * nc)
                    rhs_scaled = n_plus * nd * n_minus ** 2
                    assert lhs_scaled >= rhs_scaled
                    nonzero_bound += 1
                # Exact signed-product determinant formulas verify the common
                # d=N(kD) and the three independent norm factors.
                assert delta[(0, 1)] == scale * na * t["01"]
                assert delta[(2, 3)] == scale * na * t["23"]
                assert delta[(0, 2)] == scale * nb * t["02"]
                assert delta[(1, 3)] == scale * nb * t["13"]
                assert delta[(0, 3)] == scale * nc * t["03"]
                assert delta[(1, 2)] == scale * nc * t["12"]
                accepted += 1
    assert accepted > 0 and nonzero_bound > 0
    return accepted, nonzero_bound


if __name__ == "__main__":
    random_count, nontrivial_count, random_crosses = check_random_coordinate_fixtures()
    product_count, bound_count = check_gaussian_product_fixtures()
    varying_count, strong_count, unit_count, zero_delta_count, max_norm_count, exception_count = \
        check_varying_corrections_and_exceptions()
    relabel_count = check_exact_group_reorientations()
    print(f"PASS: {random_count} retained arbitrary-coordinate fixtures; "
          f"{nontrivial_count} have a nontrivial coprime matching factor; "
          f"{random_crosses} nonzero row-pair cross products checked.")
    print(f"PASS: {product_count} signed Gaussian-product fixtures have all six "
          "matching divisibilities, nonzero coordinates/determinants, rank-two "
          "triple rows, primitive kernel, and the general height bound.")
    print(f"PASS: {bound_count} product fixtures satisfy the exact scaled "
          "quantitative inequality when k<n_minus.")
    print(f"PASS: {varying_count} varying-correction fixtures; {strong_count} "
          f"nonvanishing quartets satisfy n_D<=8 M^2 Y^4, including {unit_count} "
          f"D-unit cases; {max_norm_count} also check the max-norm branch; universal gap checked throughout, including "
          f"{exception_count} forced vanishing exceptions.")
    print(f"PASS: {relabel_count} exact group-order/orientation/common-D "
          "normalizations.")
