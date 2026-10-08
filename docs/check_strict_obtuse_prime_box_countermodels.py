"""Exact prime-weight box countermodels, not actual short-arc examples.

Tests common-unit Gaussian realization and norm conditions for both the
minimal Walsh model and the fivefold model with distinct pair supports.
All prime certificates use trial division; the checker needs no packages.
The asymptotic dyadic prime-count argument is a separate analytic input.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt, prod


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n) + 1, 2))


def nearby_split_primes(M):
    r = 5 * (M - 1)
    X = (8 * M) ** 4  # C=1 in the literal finite fixtures.
    primes = []
    candidate = X + (1 - X) % 4
    while len(primes) < r:
        assert r * (candidate - X) < X
        if is_prime(candidate):
            primes.append(candidate)
        candidate += 4
    assert len(set(primes)) == r
    assert all(p % 4 == 1 for p in primes)
    assert r * (max(primes) - min(primes)) < min(primes)
    return X, primes


def allocation_rows(M, copies, flipped):
    assert M >= 4 and M & (M - 1) == 0
    columns = [a for _ in range(copies) for a in range(1, M)]
    rows = [[bin(x & a).count("1") % 2 for a in columns] for x in range(M)]
    if flipped:
        assert len(columns) >= M
        for x in range(M):
            rows[x][x] ^= 1
    assert all({row[j] for row in rows} == {0, 1} for j in range(len(columns)))
    return rows


def rational_rank(matrix):
    rows = [[Fraction(x) for x in row] for row in matrix]
    rank = 0
    for column in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_value = rows[rank][column]
        rows[rank] = [x / pivot_value for x in rows[rank]]
        for i in range(rank + 1, len(rows)):
            multiple = rows[i][column]
            if multiple:
                rows[i] = [x - multiple * y for x, y in zip(rows[i], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def gaussian_mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def gaussian_conj(z):
    return z[0], -z[1]


def gaussian_norm(z):
    return z[0] * z[0] + z[1] * z[1]


def gaussian_prime(p):
    for a in range(1, isqrt(p) + 1):
        b = isqrt(p - a * a)
        if b > 0 and a * a + b * b == p:
            assert gcd(a, b) == 1
            return a, b
    raise AssertionError(p)


def nearest(n, d):
    return (2 * n + d) // (2 * d)


def gaussian_gcd(a, b):
    while b != (0, 0):
        raw = gaussian_mul(a, gaussian_conj(b))
        denominator = gaussian_norm(b)
        quotient = tuple(nearest(x, denominator) for x in raw)
        multiple = gaussian_mul(quotient, b)
        a, b = b, (a[0] - multiple[0], a[1] - multiple[1])
    return a


def gaussian_rows(rows, prime_factors):
    output = []
    for row in rows:
        value = (1, 0)
        for bit, factor in zip(row, prime_factors):
            value = gaussian_mul(value, factor if bit else gaussian_conj(factor))
        output.append(value)
    return output


def check_model(M, primes, prime_factors, copies, flipped):
    rows = allocation_rows(M, copies, flipped)
    r = len(rows[0])
    primes, prime_factors = primes[:r], prime_factors[:r]
    N = prod(primes)
    points = gaussian_rows(rows, prime_factors)
    assert len(set(points)) == M
    assert all(gaussian_norm(z) == N for z in points)
    common = points[0]
    for point in points[1:]:
        common = gaussian_gcd(common, point)
    assert gaussian_norm(common) == 1
    differences = [[x - y for x, y in zip(row, rows[0])] for row in rows[1:]]
    assert rational_rank(differences) == M - 1

    if not flipped:
        signs = [[2 * bit - 1 for bit in row] for row in rows]
        assert all(sum(row[j] for row in signs) == 0 for j in range(r))
        # F^T F = M diag(log p_j), certifying the stated nonzero spectrum.
        for j in range(r):
            for k in range(r):
                assert sum(row[j] * row[k] for row in signs) == M * int(j == k)

    supports = []
    norms = []
    for i, j in combinations(range(M), 2):
        support = tuple(k for k in range(r) if rows[i][k] != rows[j][k])
        hamming = len(support)
        assert 2 * hamming - r >= 1
        P = prod(primes[k] for k in support)
        assert N % P == 0
        # Integer forms of strict obtuseness and the stronger abstract-arc
        # chord packing bound, without floating-point logarithms or sines.
        assert P * P > N
        assert P * P > 256 * (M - 1) ** 4 * N
        common_pair = gaussian_gcd(points[i], points[j])
        assert gaussian_norm(common_pair) == N // P
        chord = tuple(a - b for a, b in zip(points[i], points[j]))
        assert gaussian_norm(chord) >= 4 * (N // P)
        supports.append(support)
        norms.append(P)
    if flipped:
        assert len(set(supports)) == len(supports)
        assert len(set(norms)) == M * (M - 1) // 2
    else:
        assert len(set(supports)) == len(set(norms)) == M - 1
    return len(norms)


def main():
    total_pairs = 0
    for M in (4, 8, 16):
        X, primes = nearby_split_primes(M)
        prime_factors = [gaussian_prime(p) for p in primes]
        for p, factor in zip(primes, prime_factors):
            assert gaussian_norm(factor) == p
        first = check_model(M, primes, prime_factors, copies=1, flipped=False)
        second = check_model(M, primes, prime_factors, copies=5, flipped=True)
        total_pairs += first + second
        print(f"PASS: M={M}, support {M-1} and {5*(M-1)}, "
              f"split primes {primes[0]}..{primes[-1]} near X={X}.")
    print(f"PASS: {total_pairs} exact pair tests; affine ranks, Gram spectrum certificate,")
    print("Gaussian norms/gcds, strict norm gaps, and abstract-arc chord bounds.")
    print("PASS: every pair norm is distinct in each modified model.")
    print("Scope: the actual Gaussian arguments are not asserted to form a short arc.")


if __name__ == "__main__":
    main()
