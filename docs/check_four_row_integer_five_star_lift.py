"""Exact integer check of the four-row five-star lift and its kernel minors."""

from itertools import combinations
from math import gcd, isqrt, prod


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def add(z, w):
    return z[0] + w[0], z[1] + w[1]


def scale(a, z):
    return a * z[0], a * z[1]


def conj(z):
    return z[0], -z[1]


def norm(z):
    return z[0] ** 2 + z[1] ** 2


def product(values):
    result = (1, 0)
    for value in values:
        result = mul(result, value)
    return result


def bracket(z, w):
    return z[0] * w[1] - z[1] * w[0]


def gaussian_prime(p):
    for a in range(1, isqrt(p) + 1):
        b2 = p - a * a
        b = isqrt(b2)
        if b > 0 and b * b == b2:
            return a, b
    raise AssertionError(f"{p} is not a split prime")


def bezout(a, b):
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if old_r < 0:
        return -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def integer_kernel_basis(qs):
    """Column Euclid: return an integral basis of ker_Z[Re Q; Im Q]."""
    mat = [[z[0] for z in qs], [z[1] for z in qs]]
    change = [[int(i == j) for j in range(5)] for i in range(5)]
    for pivot in range(2):
        first = next(c for c in range(pivot, 5) if mat[pivot][c])
        for row in mat + change:
            row[pivot], row[first] = row[first], row[pivot]
        for c in range(pivot + 1, 5):
            a, b = mat[pivot][pivot], mat[pivot][c]
            if not b:
                continue
            g, s, t = bezout(a, b)
            for row in mat + change:
                x, y = row[pivot], row[c]
                row[pivot] = s * x + t * y
                row[c] = -(b // g) * x + (a // g) * y
        assert all(mat[pivot][c] == 0 for c in range(pivot + 1, 5))
    assert all(mat[r][c] == 0 for r in range(2) for c in range(2, 5))
    return [[change[i][c] for i in range(5)] for c in range(2, 5)]


def det3(rows):
    a, b, c = rows
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def determinant(rows):
    if len(rows) == 1:
        return rows[0][0]
    if len(rows) == 2:
        return rows[0][0] * rows[1][1] - rows[0][1] * rows[1][0]
    return det3(rows)


def check_correction_valuation_depth():
    """Use actual relation lattices; do not delete correction primes.

    These local fixtures test the rule for arbitrary Gaussian stars,
    including ordinary contents. They are not endpoint profiles.
    """
    p, pi = 5, (2, 1)
    barpi = conj(pi)
    units_at_p = [(1, 1), (1, 4), (2, 2), (2, 5), (3, 2)]
    assert all(norm(z) % p for z in units_at_p)
    corrections = [
        ([0, 0, 0, 0, 0], [0, 0, 0, 0, 0]),
        ([0, 1, 0, 0, 0], [0, 0, 1, 0, 0]),
        ([0, 2, 1, 3, 1], [0, 1, 2, 1, 2]),
        ([0, 0, 4, 1, 2], [0, 3, 0, 2, 1]),
    ]
    for size in (1, 2):
        for support_tuple in combinations(range(5), size):
            support = set(support_tuple)
            for exponent in (1, 4, 20):
                for kplus, kminus in corrections:
                    bplus = [exponent * (c in support) + kplus[c] for c in range(5)]
                    bminus = kminus
                    qs = [product([units_at_p[c]] + [pi] * bplus[c]
                                  + [barpi] * bminus[c]) for c in range(5)]
                    basis = integer_kernel_basis(qs)
                    for row in basis:
                        assert sum(row[c] * qs[c][0] for c in range(5)) == 0
                        assert sum(row[c] * qs[c][1] for c in range(5)) == 0
                    for r in (1, 2, 3):
                        for columns in combinations(range(5), r):
                            complement = set(range(5)) - set(columns)
                            dplus = max(0, min(bplus[c] for c in complement) - min(bplus))
                            dminus = max(0, min(bminus[c] for c in complement) - min(bminus))
                            modulus = p ** max(dplus, dminus)
                            for selected_rows in combinations(range(3), r):
                                minor = determinant([[basis[a][c] for c in columns]
                                                     for a in selected_rows])
                                assert minor % modulus == 0
                            if r <= 2 or not complement <= support:
                                assert dplus <= max(kplus)
                                assert dminus <= max(kminus)
                                assert max(dplus, dminus) <= sum(kplus) + sum(kminus)
                            else:
                                assert r == 3 and complement == support
                                assert abs(dplus - exponent) <= max(kplus)


def check_one_coordinate_reconstruction():
    """Audit Section 4 of four_row_one_coordinate_divisor_reduction.md.

    This is the existing full-support *three*-row fixture; it does not
    assert the unresolved four-row small-coordinate profile.
    """
    positive_rows = [(1469178, 1), (31686, 1), (153548, 1)]
    for sigma in (1, -1):
        rows = [(x, sigma * y) for x, y in positive_rows]
        a, signed_b = rows[0]
        b = abs(signed_b)
        n1 = norm(rows[0])
        assert isqrt(n1) == a and n1 - a * a == b * b
        assert b * b < 2 * a - 1 and gcd(a, b) == gcd(a, n1) == 1
        g = conj(rows[0])
        assert g == (a, -sigma * b)
        for row in rows[1:]:
            s = rows[0][0] * row[0] + rows[0][1] * row[1]
            delta = bracket(rows[0], row)
            z = (s, delta)
            assert mul(g, row) == z
            numerator = mul(z, conj(g))
            assert numerator[0] % n1 == numerator[1] % n1 == 0
            assert (numerator[0] // n1, numerator[1] // n1) == row
            assert (a * s - sigma * b * delta) % n1 == 0
            assert (sigma * b * s + a * delta) % n1 == 0
            assert row[1] * n1 == s * rows[0][1] + delta * rows[0][0]


def check_imaginary_content_counterexample():
    """A shared Gaussian factor survives determinant normalization, not norms."""
    h = (2, 1)
    p1 = mul(h, (-1, 2))
    p2 = mul(h, (5, 2))
    assert p1 == (-4, 3) and p2 == (8, 9)
    assert gcd(*p1) == gcd(*p2) == 1
    assert norm(p1) % 2 == norm(p2) % 2 == 1
    shared = norm(h)
    delta = bracket(p1, p2)
    t = delta // shared
    g = gcd(abs(p1[1]), abs(p2[1]))
    assert (shared, delta, t, g) == (5, -60, -12, 3)
    assert gcd(abs(p1[1]), abs(t)) == gcd(abs(p2[1]), abs(t)) == g
    virtual1 = (p1[0], p1[1] // g)
    virtual2 = (p2[0], p2[1] // g)
    assert bracket(virtual1, virtual2) == shared * (t // g)
    assert norm(virtual1) == 17 and norm(virtual2) == 73
    assert norm(virtual1) % shared and norm(virtual2) % shared


def check_actual_row_private_gcd(base_blocks):
    """Full-profile gcd/sign tests, including co-singleton overlap depth.

    These are primitive fifteen-block rows, not claimed near-axis
    endpoint fixtures. V_i itself need not remain primitive.
    """
    checks = 0
    sign_checks = 0
    collapse_checks = 0
    root_checks = 0
    content_checks = 0
    nontrivial_content = 0
    global_content_fixtures = 0
    virtual_norm_failures = 0
    for core_depth in (1, 2, 3):
        blocks = {mask: product([z] * core_depth)
                  for mask, z in base_blocks.items()}
        for correction_depth in (0, 1, 2, 4):
            corrections = {
                i: product([blocks[15 ^ (1 << (i - 1))]]
                           * correction_depth)
                for i in range(1, 5)
            }
            rows = {
                i: mul(corrections[i], product(
                    blocks[mask] for mask in range(1, 16)
                    if mask & (1 << (i - 1))))
                for i in range(1, 5)
            }
            assert all(gcd(abs(z[0]), abs(z[1])) == 1
                       and norm(z) % 2 for z in rows.values())
            reduced = {}
            for i in range(1, 5):
                for j in range(1, 5):
                    if i == j:
                        continue
                    shared = prod(norm(blocks[mask])
                                  for mask in range(1, 16)
                                  if mask & (1 << (i - 1))
                                  and mask & (1 << (j - 1)))
                    z = mul(conj(rows[i]), rows[j])
                    assert z[0] % shared == z[1] % shared == 0
                    reduced[i, j] = z[0] // shared, z[1] // shared
                    assert reduced[i, j][1] != 0
            global_content = gcd(*(abs(rows[i][1]) for i in range(1, 5)))
            assert all(gcd(global_content, norm(rows[i])) == 1
                       for i in range(1, 5))
            if global_content > 1:
                global_content_fixtures += 1
            virtual_rows = {
                i: (rows[i][0], rows[i][1] // global_content)
                for i in range(1, 5)
            }
            if global_content > 1 and any(
                norm(virtual_rows[i]) % norm(blocks[mask])
                for i in range(1, 5)
                for mask in range(1, 16)
                if mask & (1 << (i - 1))
            ):
                virtual_norm_failures += 1
            for i, j in combinations(range(1, 5), 2):
                yi, yj = rows[i][1], rows[j][1]
                tij = reduced[i, j][1]
                common = gcd(abs(yi), abs(yj))
                assert gcd(abs(yi), abs(tij)) == common
                assert gcd(abs(yj), abs(tij)) == common
                assert tij % global_content == 0
                assert gcd(abs(virtual_rows[i][1]), abs(tij // global_content)) == (
                    gcd(abs(virtual_rows[i][1]), abs(virtual_rows[j][1])))
                shared = prod(norm(blocks[mask])
                              for mask in range(1, 16)
                              if mask & (1 << (i - 1))
                              and mask & (1 << (j - 1)))
                assert bracket(virtual_rows[i], virtual_rows[j]) == (
                    shared * (tij // global_content))
                content_checks += 1
                if common > 1:
                    nontrivial_content += 1
            for i in range(1, 5):
                singleton = norm(blocks[1 << (i - 1)])
                co_singleton = norm(blocks[15 ^ (1 << (i - 1))])
                h = norm(corrections[i]) * singleton
                grouped = h * co_singleton
                quotients = []
                others = [j for j in range(1, 5) if j != i]
                for j, k in combinations(others, 2):
                    xij, tij = reduced[i, j]
                    xik, tik = reduced[i, k]
                    tjk = reduced[j, k][1]
                    numerator = tik * xij - tij * xik
                    assert numerator % tjk == 0
                    quotients.append(abs(numerator // tjk))
                assert gcd(*quotients) == grouped
                x_i, y_i = rows[i]
                assert gcd(x_i, singleton) == 1
                for candidate_x, candidate_y in (
                    (0, 0), (1, 0), (0, 1), (2, -3), (x_i, y_i)
                ):
                    outputs = [gcd(grouped, candidate_x * reduced[i, j][1]
                                   + candidate_y * reduced[i, j][0])
                               for j in others]
                    common = gcd(*outputs)
                    for j, output in zip(others, outputs):
                        assert abs(reduced[i, j][1]) % (output // common) == 0
                        collapse_checks += 1
                for j in others:
                    xij, tij = reduced[i, j]
                    assert (xij * xij + tij * tij) % grouped == 0
                    u = gcd(xij, tij)
                    reduced_modulus = grouped // gcd(grouped, u * u)
                    assert reduced_modulus * tij * tij >= grouped
                    assert gcd(tij // u, reduced_modulus) == 1
                    if reduced_modulus > 1:
                        root = (xij // u) * pow(tij // u, -1, reduced_modulus)
                        assert (root * root + 1) % reduced_modulus == 0
                    root_checks += 1
                    correct = gcd(grouped, x_i * tij + y_i * xij)
                    opposite = gcd(grouped, x_i * tij - y_i * xij)
                    assert correct == h
                    assert grouped // correct == co_singleton
                    if correct == opposite == h:
                        assert tij % singleton == 0
                    if tij % singleton != 0:
                        assert opposite != h
                        sign_checks += 1
                    checks += 1
    assert checks == 144 and sign_checks > 0
    assert collapse_checks == 720 and root_checks == 144
    assert content_checks == 72 and nontrivial_content > 0
    assert global_content_fixtures > 0 and virtual_norm_failures > 0
    return (checks, sign_checks, collapse_checks, root_checks,
            content_checks, nontrivial_content, global_content_fixtures,
            virtual_norm_failures)


def main():
    # Fifteen distinct odd split primes for the Boolean blocks; the four
    # correction primes are disjoint from them and test their exact placement.
    block_primes = [5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113, 137]
    correction_primes = [149, 157, 181, 193]
    blocks = {mask: gaussian_prime(p) for mask, p in enumerate(block_primes, 1)}
    corrections = {a: gaussian_prime(correction_primes[a - 1]) for a in range(1, 5)}
    assert len({norm(z) for z in blocks.values()}) == 15
    assert all(norm(z) % 2 == 1 and gcd(*z) == 1 for z in blocks.values())

    rows = {a: mul(corrections[a], product(blocks[mask]
            for mask in range(1, 16) if mask & (1 << (a - 1)))) for a in range(1, 5)}
    assert all(gcd(*z) == 1 for z in rows.values())

    vertices = {0: conj(blocks[15])}
    vertices.update({a: blocks[1 << (a - 1)] for a in range(1, 5)})
    edges = {(0, a): conj(blocks[15 ^ (1 << (a - 1))]) for a in range(1, 5)}
    edges.update({(a, b): blocks[(1 << (a - 1)) | (1 << (b - 1))]
                  for a, b in combinations(range(1, 5), 2)})
    qs = {0: product([vertices[0]] + [edges[0, a] for a in range(1, 5)])}
    for a in range(1, 5):
        qs[a] = product([corrections[a], vertices[a]] +
                        [edges[min(a, b), max(a, b)] for b in range(5) if b != a])

    for a in range(1, 5):
        assert mul(conj(qs[0]), qs[a]) == scale(norm(edges[0, a]), rows[a])

    residuals = {}
    minors = {}
    for a, b in combinations(range(5), 2):
        d = bracket(qs[a], qs[b])
        n = norm(edges[a, b])
        assert d % n == 0
        t = d // n
        assert t != 0
        if a == 0:
            assert t == rows[b][1]
        else:
            shared = prod(norm(blocks[mask]) for mask in range(1, 16)
                          if mask & (1 << (a - 1)) and mask & (1 << (b - 1)))
            assert bracket(rows[a], rows[b]) % shared == 0
            assert t == bracket(rows[a], rows[b]) // shared
        residuals[a, b] = t
        minors[a, b] = d

    h = gcd(*(abs(d) for d in minors.values()))
    for e, f in combinations(edges, 2):
        assert h > 0 and abs(residuals[e] * residuals[f]) % h == 0

    for a, b, c in combinations(range(5), 3):
        relation = [0] * 5
        relation[a] = minors[b, c]
        relation[b] = -minors[a, c]
        relation[c] = minors[a, b]
        content = gcd(*(abs(x) for x in relation))
        assert abs(residuals[a, b] * residuals[a, c]) % content == 0
        assert sum(relation[u] * qs[u][0] for u in range(5)) == 0
        assert sum(relation[u] * qs[u][1] for u in range(5)) == 0

    basis = integer_kernel_basis([qs[a] for a in range(5)])
    for lam in basis:
        assert sum(lam[a] * qs[a][0] for a in range(5)) == 0
        assert sum(lam[a] * qs[a][1] for a in range(5)) == 0
        for a, b in edges:
            complement = (0, 0)
            for c in range(5):
                if c not in (a, b):
                    complement = add(complement, scale(lam[c], qs[c]))
            assert all(v % norm(edges[a, b]) == 0
                       for v in mul(complement, conj(edges[a, b])))

    squares = []
    for a, b in edges:
        remaining = [c for c in range(5) if c not in (a, b)]
        minor = det3([[basis[r][c] for c in remaining] for r in range(3)])
        assert abs(minor) == abs(minors[a, b]) // h
        squares.append(minor * minor)
    assert h * h * sum(squares) == sum(d * d for d in minors.values())
    check_correction_valuation_depth()
    check_one_coordinate_reconstruction()
    check_imaginary_content_counterexample()
    (private_checks, sign_checks, collapse_checks, root_checks,
     content_checks, nontrivial_content, global_content_fixtures,
     virtual_norm_failures) = (
        check_actual_row_private_gcd(blocks))
    print("PASS: exact integer K5 lift, kernel minors with correction-depth losses, "
          "and three-row nearest-square reconstruction")
    print(f"PASS: {private_checks} actual-row private-gcd checks and "
          f"{sign_checks} sign exclusions, with core/correction valuation overlap")
    print(f"PASS: {collapse_checks} arbitrary-candidate gcd collapses and "
          f"{root_checks} valuation-safe modular roots")
    print(f"PASS: {content_checks} residual-content identities "
          f"({nontrivial_content} nontrivial pairs), "
          f"{global_content_fixtures} common-content fixtures, and "
          f"{virtual_norm_failures} virtual norm failures")


if __name__ == "__main__":
    main()
