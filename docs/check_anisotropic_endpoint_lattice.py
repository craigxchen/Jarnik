"""Exact audit for the centered anisotropic endpoint coordinates."""

from math import gcd, hypot, sqrt

from check_direct_primitive_axis_lift import family
from check_gaussian_reflection_replacement import gnorm


def bezout(a, b):
    """Return r,s with a*r+b*s=gcd(a,b), for gcd(a,b)=1."""
    old_r, r = 1, 0
    old_s, s = 0, 1
    aa, bb = a, b
    while bb:
        quotient = aa // bb
        aa, bb = bb, aa - quotient * bb
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
    if aa == -1:
        aa = 1
        old_r, old_s = -old_r, -old_s
    assert aa == 1
    return old_r, old_s


def check_endpoint(z0, rows):
    a, b = z0
    g = gcd(abs(a), abs(b))
    p, q = a // g, b // g
    A = p * p + q * q
    rho2 = gnorm(z0)
    assert rho2 == g * g * A
    h, k = bezout(p, q)
    # The Euclidean algorithm returns p*h+q*k=1.  Rotate that identity.
    r, s = -k, h
    assert p * s - q * r == 1
    B = p * r + q * s
    assert A * (r * r + s * s) - B * B == 1
    assert (B * B + 1) % A == 0

    local = []
    for z in rows:
        dx, dy = z[0] - z0[0], z[1] - z0[1]
        x = p * dx + q * dy
        y = p * dy - q * dx
        if dx == dy == 0:
            assert (x, y) == (0, 0)
            continue
        d = -x
        assert d > 0
        assert x * x + y * y == 2 * g * A * d
        assert y * y == d * (2 * g * A - d)
        assert (d + B * y) % A == 0
        assert (dx * dx + dy * dy) * A == x * x + y * y
        # Numeric values are only used to display the collapse; all identities
        # above are exact integer assertions.
        rho = sqrt(rho2)
        normal = x / sqrt(A)
        tangent_scaled = y / sqrt(A * rho)
        assert abs(tangent_scaled * tangent_scaled + 2 * normal
                   + normal * normal / rho) < 1e-12
        local.append((normal, tangent_scaled))
    return local


def main():
    maxima = []
    for index in (1, 11, 21, 31):
        U, V, P, A_block, B_block, C_block, rows = family(index)
        assert all(gnorm(z) == gnorm(rows[0]) for z in rows)
        for endpoint in rows:
            local = check_endpoint(endpoint, rows)
            maxima.append((V, max((hypot(n, t) for n, t in local), default=0.0)))
    assert maxima[-1][1] < maxima[0][1]
    print("PASS: exact centered dot/cross parabola and congruence checks.")
    print("PASS: primitive four-point Pell endpoints collapse under the "
          "anisotropic normalization.")


if __name__ == "__main__":
    main()
