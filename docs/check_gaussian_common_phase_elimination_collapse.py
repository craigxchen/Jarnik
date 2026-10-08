"""Exact checks for common-phase elimination and scalar collapse."""

from math import gcd
from random import Random


def gmul(z, w):
    return (z[0] * w[0] - z[1] * w[1],
            z[0] * w[1] + z[1] * w[0])


def gpow(z, n):
    out = (1, 0)
    for _ in range(n):
        out = gmul(out, z)
    return out


def conjugate(z):
    return (z[0], -z[1])


def composite(primes, exponents):
    out = (1, 0)
    for pi, exponent in zip(primes, exponents):
        factor = pi if exponent >= 0 else conjugate(pi)
        out = gmul(out, gpow(factor, abs(exponent)))
    return out


def allocation_checks():
    rng = Random(478321)
    primes = [(2, 1), (3, 2), (4, 1)]  # Norms 5, 13, 17.
    norms = [x * x + y * y for x, y in primes]
    units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    for _ in range(500):
        rows = 8
        allocations = [
            sorted(rng.randrange(7) for _ in range(rows))
            for _ in primes
        ]
        l = [rng.randrange(-3, 4) for _ in range(rows - 1)]
        l.append(-sum(l))
        k = [rng.randrange(-3, 4) for _ in range(rows - 1)]
        k.append(-sum(k))
        cl = [sum(l[i] * a[i] for i in range(rows)) for a in allocations]
        ck = [sum(k[i] * a[i] for i in range(rows)) for a in allocations]
        al, ak = composite(primes, cl), composite(primes, ck)
        left = gmul(al, conjugate(ak))
        diff = [x - y for x, y in zip(cl, ck)]
        residual = composite(primes, diff)
        content = 1
        for p, x, y in zip(norms, cl, ck):
            exponent = (abs(x) + abs(y) - abs(x - y)) // 2
            assert exponent >= 0
            content *= p ** exponent
        assert left == (content * residual[0], content * residual[1])
        assert left[1] == content * residual[1]
        reverse = composite(primes, [-x for x in diff])
        for ul in units:
            for uk in units:
                determinant = gmul(conjugate(gmul(ul, al)), gmul(uk, ak))[1]
                twist = gmul(conjugate(ul), uk)
                rhs = content * gmul(twist, reverse)[1]
                assert determinant == rhs


def family_checks():
    a = [2, 4, 6, 8, 10]
    period = 1
    for x in a:
        for y in a:
            if x < y:
                d = y * y - x * x
                period = period * d // gcd(period, d)
    for multiplier in (1, 3, 7):
        t = multiplier * period
        values = [(x * t, 1) for x in a]
        norms = [x * x + y * y for x, y in values]
        assert all(gcd(x, y) == 1 and (x - y) % 2 for x, y in values)
        assert all(gcd(norms[i], norms[j]) == 1
                   for i in range(len(a)) for j in range(i))
        determinants = [
            abs(values[i][0] * values[j][1]
                - values[i][1] * values[j][0])
            for i in range(len(a)) for j in range(i)
        ]
        joint = 0
        for value in determinants:
            joint = gcd(joint, value)
        assert joint == 2 * t
    print("PASS: 500 allocation collapses; primitive coprime joint-ideal family")


if __name__ == "__main__":
    allocation_checks()
    family_checks()
