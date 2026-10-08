"""Exact Gaussian and finite-sample audits of source-aligned circle maps."""

from math import atan2, cos, gcd, pi, sqrt

from check_full_cut_affine_star_extraction_obstruction import mul
from check_mobius_conductor_transfer import fixtures, least_norm


def norm(z):
    return z[0] * z[0] + z[1] * z[1]


def gaussian_valuation(z, prime):
    a, b = prime
    p = a * a + b * b
    x, y = z
    count = 0
    while True:
        num_x, num_y = x * a + y * b, y * a - x * b
        if num_x % p or num_y % p:
            return count
        x, y = num_x // p, num_y // p
        count += 1


def signed_profile(rows, prime):
    conjugate = prime[0], -prime[1]
    return [gaussian_valuation(z, prime) - gaussian_valuation(z, conjugate)
            for z in rows]


def centered_image(rows, j, scale):
    x, y = rows[j]
    images, raw_norms, contents = [], [], []
    for u, v in rows:
        dot, determinant = x * u + y * v, x * v - y * u
        X, Y = scale * dot, determinant
        content = gcd(X, Y)
        images.append((X // content, Y // content))
        raw_norms.append(X * X + Y * Y)
        contents.append(content)
    assert images[j] == (1, 0)
    assert gcd(scale * x, scale * y, -y, x) == 1
    assert scale * norm(rows[j]) != 0
    return images, raw_norms, contents


def shortest_phase_arc(rows):
    angles = sorted((2 * atan2(y, x)) % (2 * pi) for x, y in rows)
    gaps = [angles[i + 1] - angles[i] for i in range(len(angles) - 1)]
    gaps.append(angles[0] + 2 * pi - angles[-1])
    return 2 * pi - max(gaps)


def full_cut_audit():
    gaussian_primes = [(2, 1), (3, 2), (4, 1), (5, 2),
                       (6, 1), (5, 4), (7, 2)]
    prime_norms = [norm(z) for z in gaussian_primes]
    assert prime_norms == [5, 13, 17, 29, 37, 41, 53]
    cuts = [mask for mask in range(1, 15) if mask & 1]
    assert len(cuts) == 7

    # At every cut, row 0 uses pi and all other rows choose pi/bar pi
    # by cut membership. H_i is the exact relative half-angle numerator.
    rows = [(1, 0)]
    for i in range(1, 4):
        h = (1, 0)
        for cut, gaussian in zip(cuts, gaussian_primes):
            if not cut & (1 << i):
                h = mul(h, (gaussian[0], -gaussian[1]))
        rows.append(h)
    N = least_norm(rows)
    source_product = 1
    for prime, cut in zip(gaussian_primes, cuts):
        source_product *= norm(prime)
        profile = signed_profile(rows, prime)
        assert profile == [0 if cut & (1 << i) else -1 for i in range(4)]
        assert max(profile) - min(profile) == 1
    assert N == source_product == 2576450045

    # One exact source-derived map erases every old cut prime. The least
    # output radius is nevertheless much larger and has external support.
    images, raw_norms, contents = centered_image(rows, 0, 3)
    assert contents[0] == 3
    assert raw_norms[0] == 9
    N_prime = least_norm(images)
    assert N_prime == 2393873733075157
    assert gcd(N, N_prime) == 1
    for prime in gaussian_primes:
        profile = signed_profile(images, prime)
        assert max(profile) - min(profile) == 0
    return len(cuts), N, N_prime


def pell_endpoint_audit():
    rows = fixtures()[3]
    N = least_norm(rows)
    assert N == 358853785
    theta = shortest_phase_arc(rows)
    assert theta < pi / 2 and theta * N ** 0.25 < 2
    checks = 0
    for j in range(len(rows)):
        x, y = rows[j]
        f_j = norm(rows[j])
        cos_sq = []
        for u, v in rows:
            dot, determinant = x * u + y * v, x * v - y * u
            f_i = u * u + v * v
            assert dot * dot + determinant * determinant == f_j * f_i
            assert dot * dot > determinant * determinant
            cos_sq.append(dot * dot / (f_j * f_i))
        sin_sq_distance = min(cos_sq)  # Distance to the contracting line.
        assert sin_sq_distance > 0.5
        for scale in (2, 3, 5):
            images, raw_norms, contents = centered_image(rows, j, scale)
            N_prime = least_norm(images)
            theta_prime = shortest_phase_arc(images)
            assert theta_prime <= theta + 1e-12
            # Exact determinant identity and finite lattice separation.
            for i, (u, v) in enumerate(rows):
                if i == j:
                    continue
                dot, determinant = x * u + y * v, x * v - y * u
                assert raw_norms[i] == scale * scale * dot * dot + determinant * determinant
                assert raw_norms[i] >= scale * f_j * norm(rows[i])
                assert 4 * N_prime * determinant * determinant >= raw_norms[i]
                assert contents[i] > 0
                checks += 1
            # Arc-lift gap bound and the nonzero-distance derivative bound.
            assert (len(rows) - 1) / sqrt(N_prime) <= theta_prime + 1e-12
            assert scale * sin_sq_distance <= theta * sqrt(N_prime) / (len(rows) - 1)
            assert theta_prime <= theta / (scale * sin_sq_distance) + 1e-12
    return checks


def main():
    cuts, N, N_prime = full_cut_audit()
    checks = pell_endpoint_audit()
    print(f"PASS: {cuts} literal full cuts and exact transformed Gaussian widths; "
          f"all old source primes erased, N={N}, N'={N_prime} from new primes.")
    print(f"PASS: {checks} exact Pell determinant/radius comparisons and "
          "source-aligned finite-sample lift inequalities.")


if __name__ == "__main__":
    main()
