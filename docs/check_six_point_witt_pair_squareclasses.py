"""Exact Hilbert-symbol checks for the six-label pair-squareclass model."""

from itertools import combinations
from math import isqrt, prod
from random import Random


EDGES = tuple(combinations(range(6), 2))
SIGNS = (1, -1, 1, -1, 1, -1)
COMPARABLE_PRIMES = (
    1037329, 1584721, 1379513, 1501889, 1869649,
    1278713, 1700513, 1209353, 1118009, 1377881,
    1242641, 1379449, 1636121, 1207417, 1747169,
)


def is_prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p)+1))


def legendre(a, p):
    value = pow(a % p, (p-1)//2, p)
    assert value in (1, p-1)
    return 1 if value == 1 else -1


def hilbert(a, b, p):
    if p == 'infinity':
        return -1 if a < 0 and b < 0 else 1
    alpha = beta = 0
    while a % p == 0:
        a //= p
        alpha += 1
    while b % p == 0:
        b //= p
        beta += 1
    if p == 2:
        exponent = ((a-1)//2)*((b-1)//2)
        exponent += alpha*((b*b-1)//8)+beta*((a*a-1)//8)
        return (-1)**(exponent % 2)
    return ((-1)**(((p-1)//2)*alpha*beta)
            * legendre(a, p)**beta * legendre(b, p)**alpha)


def hasse(coefficients, p):
    return prod(hilbert(a, b, p) for a, b in combinations(coefficients, 2))


def check_labels(primes):
    labels = dict(zip(EDGES, primes))
    coefficients = [SIGNS[i]*prod(p for edge, p in labels.items() if i in edge)
                    for i in range(6)]
    assert prod(coefficients) == -prod(primes)**2
    local = []
    for (i, j), p in labels.items():
        crossing = prod(labels[tuple(sorted((v, k)))]
                        for v in (i, j) for k in range(6) if k not in (i, j))
        expected = legendre(-SIGNS[i]*SIGNS[j]*crossing, p)
        actual = hasse(coefficients, p)
        assert actual == expected
        local.append(actual)
    return coefficients, local


def main():
    assert all(is_prime(p) and p % 8 == 1 for p in COMPARABLE_PRIMES)
    assert max(COMPARABLE_PRIMES) < 2*min(COMPARABLE_PRIMES)
    assert all(legendre(p, q) == 1 for p, q in combinations(COMPARABLE_PRIMES, 2))
    coefficients, local = check_labels(COMPARABLE_PRIMES)
    assert local == [1]*15
    hyperbolic = SIGNS
    for place in (*COMPARABLE_PRIMES, 2, 'infinity'):
        assert hasse(coefficients, place) == hasse(hyperbolic, place)

    # Test the sign-sensitive formula without the p=1 mod 8 restriction.
    pool = [p for p in range(3, 300) if is_prime(p)]
    random = Random(619)
    failures = 0
    for _ in range(40):
        _, local = check_labels(random.sample(pool, 15))
        failures += local.count(-1)
    assert failures > 0
    print('PASS: 15 comparable prime labels, 105 mutual residue symbols, '
          'all local hyperbolic invariants, and 600 general edge formulas.')


if __name__ == '__main__':
    main()
