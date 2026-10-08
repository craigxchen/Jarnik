"""Exact certificate for the finite three-row interpolation fixture."""

from fractions import Fraction
from itertools import combinations
from math import gcd


def main():
    g, r = 17, 4
    nodes = (1385, 1482, 1747)
    a, b, c = nodes
    gaps = (b - a, c - b, c - a)
    assert gaps == (97, 265, 362)
    s = (r*r + 1) // g
    assert (r*r + 1) % g == 0 and s == 1
    x = tuple(r + g*node for node in nodes)
    norm = tuple(value*value + 1 for value in x)
    f = tuple(value // g for value in norm)
    assert all(value % g == 0 for value in norm)
    A, B, C = gaps
    u, v, w = f[0] // (A*C), f[1] // (A*B), f[2] // (B*C)
    assert (f[0] % (A*C), f[1] % (A*B), f[2] % (B*C)) == (0, 0, 0)
    assert (u, v, w) == (929, 1453, 541)
    assert g == u-v+w
    rp = r + g*a
    assert 2*rp == w*B-u*A-g*C
    assert rp*rp + 1 == g*u*A*C
    assert Fraction(2*(u*B+w*A), C)-v == Fraction(35669, 181)
    assert norm == (554555402, 634939205, 882268210)
    assert (gcd(norm[0], norm[1]), gcd(norm[1], norm[2]),
            gcd(norm[0], norm[2])) == (g*A, g*B, g*C)

    factors = (g, A, B, C, u, v, w)
    assert all(gcd(left, right) == 1 for left, right in combinations(factors, 2))
    squares = ((4, 1), (9, 4), (16, 3), (19, 1),
               (23, 20), (38, 3), (21, 10))
    assert all(value == p*p+q*q and gcd(p, q) == 1
               for value, (p, q) in zip(factors, squares))

    # A second fixture checks Gaussian orientations rather than norms alone.
    g2, r2, nodes2 = 109, -33, (134, 291, 304)
    x2 = tuple(g2*node+r2 for node in nodes2)
    norm2 = tuple(value*value+1 for value in x2)
    A2, B2, C2 = (nodes2[1]-nodes2[0],
                  nodes2[2]-nodes2[1], nodes2[2]-nodes2[0])
    f2 = tuple(value//g2 for value in norm2)
    assert all(value % g2 == 0 for value in norm2)
    u2, v2, w2 = (f2[0]//(A2*C2), f2[1]//(A2*B2),
                  f2[2]//(B2*C2))
    assert (A2, B2, C2) == (157, 13, 170)
    assert (u2, v2, w2) == (73, 4513, 4549)
    assert (x2, norm2) == ((14573, 31686, 33103),
                           (212372330, 1004002597, 1095808610))
    assert g2 == u2-v2+w2
    assert all(f2[i] == (10*nodes2[i]-3)**2+(3*nodes2[i]-1)**2
               for i in range(3))
    odd_factors = (g2, A2, B2, C2//2, u2, v2, w2)
    assert all(gcd(left, right) == 1
               for left, right in combinations(odd_factors, 2))
    odd_squares = ((10, 3), (11, 6), (3, 2), (9, 2),
                   (8, 3), (48, 47), (65, 18))
    assert all(value == p*p+q*q and gcd(p, q) == 1
               for value, (p, q) in zip(odd_factors, odd_squares))

    def mul(z, w):
        return z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0]

    def divides(z, w):
        den = w[0]*w[0]+w[1]*w[1]
        return ((z[0]*w[0]+z[1]*w[1]) % den == 0
                and (z[1]*w[0]-z[0]*w[1]) % den == 0)

    common = (10, 3)
    assert mul(common, (-3, 1)) == (r2, 1)
    assert all(divides((value, 1), common) for value in x2)
    for (i, j), gap, private in zip(((0, 1), (1, 2), (0, 2)),
                                     (A2, B2, C2),
                                     ((-11, 6), (3, -2), (-11, -7))):
        candidate = mul(common, private)
        assert candidate[0]**2+candidate[1]**2 == g2*gap
        assert divides((x2[i], 1), candidate)
        assert divides((x2[j], 1), candidate)
        assert gcd(norm2[i], norm2[j]) == g2*gap
    assert gcd(gcd(*norm2[:2]), norm2[2]) == g2

    for t in range(2, 101):
        gt = 2*t*t+2*t+1
        rt = -2*t-1
        assert rt*rt+1 == 2*gt
        assert 2*abs(rt) <= gt
        assert gt == (t+1)**2 + t*t
    print("PASS: two three-row fixtures, Gaussian gcd orientations, Vieta failure, and centered-root family.")


if __name__ == "__main__":
    main()
