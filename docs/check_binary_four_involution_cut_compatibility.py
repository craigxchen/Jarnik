"""Exact local cut-compatibility checks for four-point involutions."""

from math import gcd

from check_gaussian_reflection_replacement import gconj, gmul, gnorm, ggcd, gdivexact


def valuation(z, pi):
    a, b = z
    c, d = pi
    p = gnorm(pi)
    out = 0
    while a or b:
        x, y = a * c + b * d, b * c - a * d
        if x % p or y % p:
            break
        a, b = x // p, y // p
        out += 1
    return out


def primitive(z):
    return gdivexact(z, (gcd(z[0], z[1]), 0))


def power(z, exponent):
    out = (1, 0)
    for _ in range(exponent):
        out = gmul(out, z)
    return out


def phase_row(z, anchor, radius_squared):
    h = gmul(z, gconj(anchor))
    h = (h[0] + radius_squared, h[1])
    if h == (0, 0):
        return (1, 0)
    return primitive(gdivexact(h, ggcd(h, gconj(h))))


def involution(rows, permutation):
    equations = []
    for h, target_index in zip(rows, permutation):
        x, y = h
        u, v = rows[target_index]
        # For [[a,b],[c,-a]], collinearity with (u,v) is
        # a(xv+yu) + b(yv) - c(xu) = 0.
        equations.append((x * v + y * u, y * v, -x * u))

    for i in range(4):
        for j in range(i + 1, 4):
            r, s = equations[i], equations[j]
            a = r[1] * s[2] - r[2] * s[1]
            b = r[2] * s[0] - r[0] * s[2]
            c = r[0] * s[1] - r[1] * s[0]
            if a or b or c:
                break
        if a or b or c:
            break
    g = gcd(a, b, c)
    matrix = (a // g, b // g, c // g, -a // g)
    for h, target_index in zip(rows, permutation):
        x, y = h
        u, v = rows[target_index]
        aa, bb, cc, dd = matrix
        assert (aa * x + bb * y) * v == (cc * x + dd * y) * u
    assert matrix[0] + matrix[3] == 0
    return matrix


def coefficient_interval(matrix, pi):
    a, b, c, d = matrix
    # U=2 alpha and V=2 beta for the real matrix action.
    u, v = (a + d, c - b), (a - d, c + b)
    assert u != (0, 0) and v != (0, 0)
    endpoints = (
        valuation(v, pi) - valuation(u, pi),
        valuation(gconj(u), pi) - valuation(gconj(v), pi),
    )
    return tuple(sorted(endpoints))


def clipped_width(source_interval, coefficient_interval):
    lo = max(source_interval[0], coefficient_interval[0])
    hi = min(source_interval[1], coefficient_interval[1])
    return max(0, hi - lo)


def fixture(pattern):
    pi, rho, sigma = (2, 1), (3, 2), (4, 1)
    assignments = ((0, 1, 0, 1), (0, 0, 1, 1))
    rows = []
    for at_p, at_rho, at_sigma in zip(pattern, *assignments):
        z = power(pi if at_p else gconj(pi), 2)
        z = gmul(gmul(z, rho if at_rho else gconj(rho)),
                 sigma if at_sigma else gconj(sigma))
        rows.append(z)
    radius_squared = gnorm(rows[0])
    half_angles = [phase_row(z, rows[0], radius_squared) for z in rows]
    return rows, half_angles, radius_squared


def audit(pattern, expected_source, expected_widths, expected_targets):
    pi = (2, 1)
    rows, half_angles, radius_squared = fixture(pattern)
    assert all(gnorm(z) == radius_squared for z in rows)
    signed = [valuation(h, pi) - valuation(h, gconj(pi)) for h in half_angles]
    assert signed == expected_source
    source_interval = (min(signed), max(signed))
    permutations = ((1, 0, 3, 2), (2, 3, 0, 1), (3, 2, 1, 0))
    widths = []
    targets = []
    for permutation in permutations:
        matrix = involution(half_angles, permutation)
        interval = coefficient_interval(matrix, pi)
        widths.append(clipped_width(source_interval, interval))
        targets.append(tuple(signed[j] - signed[permutation[0]] for j in permutation))
    assert widths == expected_widths
    assert targets == expected_targets
    return widths


def main():
    three_one = audit(
        (0, 0, 0, 1),
        [0, 0, 0, 2],
        [0, 0, 0],
        [(0, 0, 2, 0), (0, 2, 0, 0), (0, -2, -2, -2)],
    )
    two_two = audit(
        (0, 0, 1, 1),
        [0, 0, 2, 2],
        [2, 2, 2],
        [(0, 0, 2, 2), (0, 0, -2, -2), (0, 0, -2, -2)],
    )
    print(f"PASS: 3|1 widths {three_one}; 2|2 widths {two_two}; "
          "all exact trace-zero and valuation checks passed.")


if __name__ == "__main__":
    main()
