"""Exact audit of the Gaussian reflection replacement identity."""

from itertools import combinations
from math import gcd, isqrt


def gadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def gsub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1],
            a[0] * b[1] + a[1] * b[0])


def gconj(a):
    return (a[0], -a[1])


def gnorm(a):
    return a[0] * a[0] + a[1] * a[1]


def gdivexact(a, b):
    den = gnorm(b)
    x = (a[0] * b[0] + a[1] * b[1])
    y = (a[1] * b[0] - a[0] * b[1])
    assert x % den == 0 and y % den == 0
    return (x // den, y // den)


def gdivrem(a, b):
    den = gnorm(b)
    x = a[0] * b[0] + a[1] * b[1]
    y = a[1] * b[0] - a[0] * b[1]
    # Nearest integer quotients, with ties rounded upward, using integers only.
    q = ((2*x + den)//(2*den), (2*y + den)//(2*den))
    return gsub(a, gmul(q, b))


def ggcd(a, b):
    while b != (0, 0):
        a, b = b, gdivrem(a, b)
    return a


def gdivides(a, b):
    return gdivexact(b, a)


def normalize(g):
    if g == (0, 0):
        return g
    # Any associate is sufficient for divisibility and norm assertions.
    return g


def gcd_many(values):
    out = (0, 0)
    for value in values:
        out = ggcd(out, value)
    return normalize(out)


def points_on_circle(n):
    bound = isqrt(n)
    return [(x, y) for x in range(-bound, bound + 1)
            for y in range(-bound, bound + 1) if x * x + y * y == n]


def check(zs, l, i, j):
    n = gnorm(zs[0])
    assert all(gnorm(z) == n for z in zs)
    assert gcd_many(zs) in ((1, 0), (-1, 0), (0, 1), (0, -1))
    retained = [z for k, z in enumerate(zs) if k != l]
    G = gcd_many(retained)
    h = ggcd(zs[l], gmul(zs[i], zs[j]))
    d = gdivexact(zs[l], h)
    A = gdivexact(gmul(zs[i], zs[j]), h)
    assert gcd_many((d, A)) in ((1, 0), (-1, 0), (0, 1), (0, -1))
    # The useful divisibility is G^2 | A; check it directly.
    gdivexact(A, gmul(G, G))
    assert gdivides(G, A)
    assert gnorm(gcd_many([gmul(d, G), A])) == gnorm(G)

    znew = [gdivexact(gmul(d, z), G) if k != l else gdivexact(A, G)
            for k, z in enumerate(zs)]
    assert gcd_many(znew) in ((1, 0), (-1, 0), (0, 1), (0, -1))
    nnew = gnorm(znew[0])
    assert all(gnorm(z) == nnew for z in znew)
    assert nnew * gnorm(G) == n * gnorm(d)

    # The reverse quotient is d^2*h/G; clear its denominator symbolically.
    reverse_numerator = gmul(gmul(d, d), h)
    cleared = [gmul(G, z) if k != l else reverse_numerator
               for k, z in enumerate(znew)]
    assert gnorm(gcd_many(cleared)) == gnorm(d)
    assert [gdivexact(z, d) for z in cleared] == list(zs)


def main():
    cases = []
    for n in range(1, 101):
        pts = points_on_circle(n)
        if len(pts) < 3:
            continue
        for zs in combinations(pts, min(4, len(pts))):
            if gcd_many(zs) not in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                continue
            for l in range(len(zs)):
                retained = [k for k in range(len(zs)) if k != l]
                i = retained[0]
                for j in retained[:2]:
                    check(zs, l, i, j)
                    cases.append(1)
    print(f"PASS: {len(cases)} primitive equal-norm reflection cases")


if __name__ == "__main__":
    main()
