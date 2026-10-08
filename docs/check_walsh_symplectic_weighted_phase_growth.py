"""Exact symplectic incidences, physical characters, and weighted heights."""

from collections import Counter
from fractions import Fraction
from math import prod

from check_strict_obtuse_prime_box_countermodels import (
    gaussian_prime, gaussian_mul, gaussian_conj, gaussian_norm, gaussian_gcd,
)


def parity(x):
    return bin(x).count('1') % 2


def symplectic_j(x, t):
    y = 0
    for k in range(0, t, 2):
        y |= ((x >> k) & 1) << (k+1)
        y |= ((x >> (k+1)) & 1) << k
    return y


def bilinear(x, y, t):
    return parity(x & symplectic_j(y, t))


def lagrangians(t):
    spaces = {frozenset({0})}
    for _ in range(t//2):
        expanded = set()
        for space in spaces:
            for u in range(1, 1 << t):
                if u not in space and all(bilinear(u, v, t) == 0 for v in space):
                    expanded.add(frozenset(space | {u ^ v for v in space}))
        spaces = expanded
    return sorted(spaces, key=lambda space: tuple(sorted(space)))


def certificates(t):
    M, h = 1 << t, 1 << (t//2)
    for V in lagrangians(t):
        assert len(V) == h and all(not bilinear(u, v, t) for u in V for v in V)
        seen = set(V)
        for x0 in range(1, M):
            if x0 in seen:
                continue
            P = {x0 ^ v for v in V}
            seen.update(P)
            c = {x: (-1)**bilinear(x0, x ^ x0, t) for x in P}
            assert 0 not in P and len(c) == h and sum(c.values()) == 0
            yield V, x0, c


def assignments(t):
    M = 1 << t
    # Nonzero rows use the first physical copy; row zero uses a second copy.
    result = {x: symplectic_j(x, t)-1 for x in range(1, M)}
    result[0] = M-1
    assert len(set(result.values())) == M
    return result


def physical_character(t, copies, c):
    M, h = 1 << t, 1 << (t//2)
    old = [sum(value*((-1)**parity(a & x)) for x, value in c.items())//2
           for a in range(1, M)]
    assert all(abs(value) in (0, h//2) for value in old)
    inverse = {j: x for x, j in assignments(t).items()}
    values, reductions = [], 0
    for j in range(copies*(M-1)):
        a = 1+j % (M-1)
        value = old[a-1]
        if j in inverse:
            x = inverse[j]
            after = value-c.get(x, 0)*((-1)**parity(a & x))
            if x in c:
                assert abs(value)-abs(after) == 1
                reductions += 1
            else:
                assert value == after
            value = after
        values.append(value)
    assert reductions == h and sum(abs(value) for value in values) == copies*M//2-h
    return old, values


def incidence_and_weights(t):
    M, h, copies = 1 << t, 1 << (t//2), 5
    rows, labels = Counter(), Counter()
    total_coefficients = [0]*(copies*(M-1))
    count = reductions = 0
    for V, x0, c in certificates(t):
        old, physical = physical_character(t, copies, c)
        selected = {a for a, value in enumerate(old, 1) if value}
        assert selected == {symplectic_j(x, t) for x in c}
        assert all(all(parity(a & v) == bilinear(x0, v, t) for v in V)
                   for a in selected)
        rows.update(c.keys())
        labels.update(selected)
        for j, value in enumerate(physical):
            total_coefficients[j] += abs(value)
        count += 1
        reductions += h
    expected = Fraction(count*h, M-1)
    assert expected.denominator == 1
    assert set(rows) == set(labels) == set(range(1, M))
    assert set(rows.values()) == set(labels.values()) == {expected.numerator}
    flipped = {j for x, j in assignments(t).items() if x != 0}
    for j, value in enumerate(total_coefficients):
        predicted = expected*h/2-expected*int(j in flipped)
        assert value == predicted
    # These arbitrary positive integers stand for log weights; their exact
    # coefficient identity therefore applies to arbitrary real prime logs.
    for shift in (1, 7, 53):
        weights = [shift+(17*j*j+3*j) % 101 for j in range(len(total_coefficients))]
        W, F = sum(weights), sum(weights[j] for j in flipped)
        actual = Fraction(sum(a*w for a, w in zip(total_coefficients, weights)), count)
        assert actual-Fraction(W, 2) == Fraction(W-2*h*F, 2*(M-1))
    return count, reductions, expected.numerator


def gaussian_product(values):
    result = (1, 0)
    for value in values:
        result = gaussian_mul(result, value)
    return result


def literal_gaussian_fixture():
    t, copies, M = 2, 5, 4
    primes = [5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113, 137]
    assert len(primes) == copies*(M-1)
    factors = [gaussian_prime(p) for p in primes]
    assert all(gaussian_norm(z) == p for z, p in zip(factors, primes))
    units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    common = (2, 1)
    assigned = assignments(t)
    points = []
    for x in range(M):
        factors_x = []
        for j, factor in enumerate(factors):
            sign = (-1)**parity(x & (1+j % (M-1)))
            if j == assigned[x]:
                sign *= -1
            factors_x.append(factor if sign == 1 else gaussian_conj(factor))
        points.append(gaussian_mul(gaussian_mul(common, units[x]), gaussian_product(factors_x)))
    N0 = prod(primes)
    assert all(gaussian_norm(z) == 5*N0 for z in points)
    shared = points[0]
    for z in points[1:]:
        shared = gaussian_gcd(shared, z)
    assert gaussian_norm(shared) == 5
    norms, count = [], 0
    for V, x0, c in certificates(t):
        old, physical = physical_character(t, copies, c)
        beta = gaussian_product(
            factor if exponent > 0 else gaussian_conj(factor)
            for factor, exponent in zip(factors, physical)
            for _ in range(abs(exponent))
        )
        expected = prod(p**abs(exponent) for p, exponent in zip(primes, physical))
        assert gaussian_norm(beta) == expected > 1
        assert gaussian_norm(gaussian_gcd(beta, gaussian_conj(beta))) == 1
        assert beta[0] and beta[1] and abs(beta[0]) != abs(beta[1])
        positive = gaussian_product(points[x] for x, sign in c.items() if sign == 1)
        negative = gaussian_product(points[x] for x, sign in c.items() if sign == -1)
        unit = units[sum(x*sign for x, sign in c.items()) % 4]
        assert gaussian_mul(positive, gaussian_conj(beta)) == gaussian_mul(
            gaussian_mul(unit, beta), negative)
        norms.append(expected)
        count += 1
    F_product = prod(primes[j] for x, j in assigned.items() if x != 0)
    assert prod(norms)*F_product**2 == N0**2
    return count


def main():
    total = reductions = 0
    for t, expected in ((4, (45, 180, 12)), (6, (945, 7560, 120))):
        result = incidence_and_weights(t)
        assert result == expected
        total += result[0]
        reductions += result[1]
    gaussian = literal_gaussian_fixture()
    print(f'PASS: {total} exact Lagrangian-coset certificates; {reductions} selected '
          f'physical reductions; six arbitrary-weight identities; {gaussian} literal '
          'Gaussian signed products and exact multiplicative average, with row units and common content.')


if __name__ == '__main__':
    main()
