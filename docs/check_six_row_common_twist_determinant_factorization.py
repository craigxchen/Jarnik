"""Exact six-row common-twist factorization and full-cut coefficient audit."""

from itertools import combinations
from math import gcd
from random import Random


PRIMES = ((2, 1, 5), (3, 2, 13), (4, 1, 17))
CHARACTERS = tuple(
    tuple(1 if i in chosen else -1 for i in range(6))
    for chosen in combinations(range(6), 3) if 0 in chosen
)


def mul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def scale(n, z):
    return (n*z[0], n*z[1])


def conjugate(z):
    return (z[0], -z[1])


def norm(z):
    return z[0]*z[0]+z[1]*z[1]


def primitive(exponents):
    result = (1, 0)
    for (real, imag, _), exponent in zip(PRIMES, exponents):
        factor = (real, imag if exponent >= 0 else -imag)
        for _ in range(abs(exponent)):
            result = mul(result, factor)
    return result


def signed_exponents(character, allocations):
    return tuple(sum(character[i]*a[i] for i in range(6))
                 for a in allocations)


def half_height(exponents, parity):
    return int_product(p**((abs(c)-e)//2)
                       for (_, _, p), c, e in zip(PRIMES, exponents, parity))


def int_product(items):
    result = 1
    for item in items:
        result *= item
    return result


def check_allocations(allocations):
    parity = tuple(sum(a) % 2 for a in allocations)
    twist = int_product(p for (_, _, p), e in zip(PRIMES, parity) if e)
    data = []
    for character in CHARACTERS:
        c = signed_exponents(character, allocations)
        assert all(x % 2 == e for x, e in zip(c, parity))
        A = primitive(c)
        h = half_height(c, parity)
        assert norm(A) == twist*h*h
        data.append((c, A, h))
    count = 0
    for (c, A, h), (d, B, k) in combinations(data, 2):
        g = gcd(h, k)
        plus = tuple((x+y)//2 for x, y in zip(c, d))
        minus = tuple((x-y)//2 for x, y in zip(c, d))
        assert all((x+y) % 2 == (x-y) % 2 == 0 for x, y in zip(c, d))
        P, M = primitive(plus), primitive(minus)
        assert mul(P, M) == scale(k//g, A)
        assert mul(P, conjugate(M)) == scale(h//g, B)
        assert k*A[0]-h*B[0] == -2*g*P[1]*M[1]
        assert k*A[0]+h*B[0] == 2*g*P[0]*M[0]
        assert norm(P)*norm(M) == twist*(h*k//g)**2
        for x, y, u, v in zip(c, d, plus, minus):
            assert abs(u)+abs(v) == max(abs(x), abs(y))
        count += 1
    assert count == 45


def check_full_cuts():
    cuts = tuple(chosen for r in (1, 2, 3)
                 for chosen in combinations(range(6), r)
                 if r < 3 or 0 in chosen)
    assert len(cuts) == 31
    height = lambda L: sum(abs(sum(L[i] for i in cut)) for cut in cuts)
    for lam, mu in combinations(CHARACTERS, 2):
        plus = tuple((x+y)//2 for x, y in zip(lam, mu))
        minus = tuple((x-y)//2 for x, y in zip(lam, mu))
        assert sorted((sum(map(abs, plus)), sum(map(abs, minus)))) == [2, 4]
        h_short, h_long = sorted((height(plus), height(minus)))
        assert (h_short, h_long) == (16, 24)
        max_height = sum(max(abs(sum(lam[i] for i in cut)),
                             abs(sum(mu[i] for i in cut))) for cut in cuts)
        assert max_height == 40


def main():
    rng = Random(20260921)
    for _ in range(40):
        allocations = []
        for _prime in PRIMES:
            exponent = rng.randrange(1, 6)
            allocations.append(tuple(rng.randrange(exponent+1) for _ in range(6)))
        check_allocations(allocations)
    check_full_cuts()
    print("PASS: 40 nested-allocation fixtures, 1,800 exact Gaussian pair factorizations")
    print("PASS: all 45 pairs on the 31-cut profile have scalar heights 16,24 and max height 40")


if __name__ == "__main__":
    main()
