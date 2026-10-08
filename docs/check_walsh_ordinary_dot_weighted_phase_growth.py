"""Exact parity certificates and restricted canonical Gaussian sources."""

from fractions import Fraction
from math import isqrt, prod

from check_walsh_symplectic_weighted_phase_growth import (
    parity, gaussian_product,
)
from check_strict_obtuse_prime_box_countermodels import (
    gaussian_prime, gaussian_mul, gaussian_conj, gaussian_norm, gaussian_gcd,
)


def even_cosets(t):
    M = 1 << t
    even = [x for x in range(M) if parity(x) == 0]
    spaces = {frozenset({0})}
    for _ in range(t // 2):
        expanded = set()
        for V in spaces:
            for u in even:
                if u not in V and all(parity(u & v) == 0 for v in V):
                    expanded.add(frozenset(V | {u ^ v for v in V}))
        spaces = expanded
    for V in spaces:
        seen = set(V)
        for x0 in even:
            if x0 in seen:
                continue
            P = {x0 ^ v for v in V}
            seen.update(P)
            c = {x: (-1)**parity(x0 & (x ^ x0)) for x in P}
            assert 0 not in c and sum(c.values()) == 0
            yield c


def odd_pairs(t):
    assert t % 2 == 0
    w = (1 << t) - 1
    for x in range(1 << t):
        if parity(x) and x < (x ^ w):
            yield {x: 1, x ^ w: -1}


def physical_columns(t, copies, canonical=False):
    M = 1 << t
    flips = ({x: x for x in range(M)} if canonical else
             {**{x: x-1 for x in range(1, M)}, 0: M-1})
    return [(a, {x for x, j in flips.items() if j == copy*(M-1)+a-1})
            for copy in range(copies) for a in range(1, M)]


def character(c, columns):
    values = []
    for a, flips in columns:
        raw = sum(value*((-1)**(parity(a & x)+int(x in flips)))
                  for x, value in c.items())
        assert raw % 2 == 0
        values.append(raw//2)
    return values


def verify_weights(t):
    M, h, w = 1 << t, 1 << (t//2), (1 << t)-1
    columns = physical_columns(t, 5)
    certs = list(even_cosets(t))
    totals = [0]*len(columns)
    for c in certs:
        assert len(c) == h
        for j, v in enumerate(character(c, columns)):
            totals[j] += abs(v)
    even_star = set(x for x in range(1, M) if parity(x) == 0 and x != w)
    for j, (a, flips) in enumerate(columns):
        if t % 2 == 0:
            q = Fraction(2*h, M-4)
            expected = q*(Fraction(h, 2)*int(a in even_star)
                          -len(flips & even_star))
        else:
            even_nonzero = set(x for x in range(1, M) if parity(x) == 0)
            q = Fraction(h, M//2-1)
            expected = q*(Fraction(h, 2)*int(a != w)
                          -len(flips & even_nonzero))
        assert Fraction(totals[j], len(certs)) == expected
    pairs = list(odd_pairs(t)) if t % 2 == 0 else []
    for c in pairs:
        values = character(c, columns)
        for v, (a, flips) in zip(values, columns):
            assert abs(v) == int(parity(a) == 1)-len(flips & c.keys())
    # Exact coefficient identities imply each arbitrary-real-weight formula.
    for shift in (1, 37):
        weights = [shift+(13*j*j+17*j)%101 for j in range(len(columns))]
        avg = Fraction(sum(v*wgt for v, wgt in zip(totals, weights)), len(certs))
        if t % 2 == 0:
            A = sum(wgt for wgt, (a, _) in zip(weights, columns) if a in even_star)
            Fe = sum(wgt*len(flips & even_star)
                     for wgt, (_, flips) in zip(weights, columns))
            assert avg == Fraction(M*A-2*h*Fe, M-4)
        else:
            active = sum(wgt for wgt, (a, _) in zip(weights, columns) if a != w)
            Fe = sum(wgt*len(flips & even_nonzero)
                     for wgt, (_, flips) in zip(weights, columns))
            assert avg == Fraction((M//2)*active-2*h*Fe, 2*(M//2-1))
    return len(certs), len(pairs)


def split_primes(n):
    result, p = [], 5
    while len(result) < n:
        if all(p % d for d in range(2, isqrt(p)+1)):
            result.append(p)
        p += 4
    return result


def literal(t, canonical=False):
    columns = physical_columns(t, 5, canonical)
    primes = split_primes(len(columns))
    factors = [gaussian_prime(p) for p in primes]
    # Mix Gaussian orientations; norms and incidence weights remain fixed.
    factors = [gaussian_conj(z) if j % 3 == 1 else z for j, z in enumerate(factors)]
    M, common = 1 << t, (2, 1)
    units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    points = []
    for x in range(M):
        row = gaussian_product(z if (parity(a & x)+int(x in flips)) % 2 == 0
                               else gaussian_conj(z)
                               for z, (a, flips) in zip(factors, columns))
        points.append(gaussian_mul(gaussian_mul(common, units[x % 4]), row))
    assert all(gaussian_norm(z) == 5*prod(primes) for z in points)
    reduced_t = t-1 if canonical else t
    certs = list(even_cosets(reduced_t))
    if reduced_t % 2 == 0:
        certs += list(odd_pairs(reduced_t))
    if canonical:
        certs = [{2*x: value for x, value in c.items()} for c in certs]
        # Every certificate avoids zero; the exceptional constant-class flip
        # at zero therefore disappears even if the class has huge prime weight.
        assert all(0 not in c for c in certs)
    for c in certs:
        values = character(c, columns)
        beta = gaussian_product(z if v > 0 else gaussian_conj(z)
                                for z, v in zip(factors, values)
                                for _ in range(abs(v)))
        assert gaussian_norm(beta) == prod(p**abs(v) for p, v in zip(primes, values)) > 1
        assert gaussian_norm(gaussian_gcd(beta, gaussian_conj(beta))) == 1
        assert beta[0] and beta[1] and abs(beta[0]) != abs(beta[1])
        positive = gaussian_product(points[x] for x, value in c.items() if value == 1)
        negative = gaussian_product(points[x] for x, value in c.items() if value == -1)
        unit = units[sum(x*value for x, value in c.items()) % 4]
        assert gaussian_mul(positive, gaussian_conj(beta)) == gaussian_mul(
            gaussian_mul(unit, beta), negative)
        if canonical:
            assert all(v == 0 for v, (a, _) in zip(values, columns) if a == 1)
    return len(certs)


def main():
    cosets = pairs = 0
    for t in (3, 4, 5, 6):
        a, b = verify_weights(t)
        cosets += a
        pairs += b
    count = sum(literal(t) for t in (3, 4))
    count += sum(literal(t, canonical=True) for t in (4, 5))
    print(f'PASS: {cosets} even-coset and {pairs} odd-pair certificates; '
          'eight arbitrary-weight identities; '
          f'{count} literal Gaussian signed products, including odd/even dimensions, '
          'mixed orientations, row units, common content, and canonical restrictions.')


if __name__ == '__main__':
    main()
