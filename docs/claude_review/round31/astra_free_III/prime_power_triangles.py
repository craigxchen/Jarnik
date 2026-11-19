"""Exact triangle test for parents with prime exponents 0, 1, or 2."""
from collections import defaultdict
from itertools import combinations, product
from math import gcd, isqrt, prod

from smooth_divisor_census import phase_compatible


def multipliers(primes, capacity, index, bound):
    if index == len(primes):
        yield 1
        return
    factor = 1
    for _ in range(capacity[index] + 1):
        if factor > bound:
            break
        for tail in multipliers(primes, capacity, index + 1, bound // factor):
            yield factor * tail
        factor *= primes[index]


def check(primes):
    parents = defaultdict(list)
    for exponents in product(range(3), repeat=len(primes)):
        d = prod(p**e for p, e in zip(primes, exponents))
        s = isqrt(d)
        offset = d - s*s
        t = isqrt(offset)
        if not t or t*t != offset:
            continue
        bound = d // (256 * offset**2)
        if bound:
            capacity = [2-e for e in exponents]
            for m in multipliers(primes, capacity, 0, bound):
                parents[d*m].append((d, exponents, gcd(s, t)))
    total = imprimitive = repeated_prime = 0
    different = set()
    for parent, candidates in parents.items():
        for triple in combinations(candidates, 3):
            if not all(sum(local) == 2*max(local)
                       for local in zip(*(exps for _, exps, _ in triple))):
                continue
            ds = tuple(sorted(d for d, _, _ in triple))
            assert phase_compatible(*ds)
            different.add(ds)
            total += 1
            imprimitive += any(g > 1 for _, _, g in triple)
            repeated_prime += any(parent % (p*p) == 0 for p in primes)
    print('Parents tested:', 3**len(primes))
    print('Profile-compatible parent/triple incidences:', total)
    print('Different divisor triples:', len(different))
    print('Incidences at nonsquarefree parents:', repeated_prime)
    print('Incidences with an imprimitive factor:', imprimitive)


if __name__ == '__main__':
    check([5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 109, 157, 229])
