"""Exact arithmetic and phase bookkeeping for the full-prime-rank lemma.

This does not numerically test, or replace, Matveev's or Roth's theorem.
Run with the Python standard library only.
"""

from itertools import combinations, product
from math import atan2, factorial, gcd, isqrt, log, pi
from random import Random


def det(matrix):
    if not matrix:
        return 1
    return sum((-1) ** j * value * det([
        row[:j] + row[j + 1:] for row in matrix[1:]
    ]) for j, value in enumerate(matrix[0]))


def adj(matrix):
    r = len(matrix)
    return [[(-1) ** (i + j) * det([
        row[:i] + row[i + 1:]
        for k, row in enumerate(matrix) if k != j
    ]) for j in range(r)] for i in range(r)]


def matmul(a, b):
    return [[sum(x * y for x, y in zip(row, col))
             for col in zip(*b)] for row in a]


def normalize(rows):
    minima = list(map(min, zip(*rows)))
    out = [[value - low for value, low in zip(row, minima)] for row in rows]
    widths = list(map(max, zip(*out)))
    return out, widths


def audit_matrix(rows):
    rows, widths = normalize(rows)
    r = len(widths)
    matrix = [[value - base for value, base in zip(row, rows[0])]
              for row in rows[1:]]
    delta = det(matrix)
    if delta == 0:
        return False
    h, ell = max(widths), factorial(r)
    jmat = adj(matrix)
    assert matmul(jmat, matrix) == [
        [delta * (i == j) for j in range(r)] for i in range(r)]
    assert abs(delta) <= ell * h ** r
    assert all(sum(map(abs, row)) <= ell * h ** (r - 1) for row in jmat)
    centered = [[2 * value - width for value, width in zip(row, widths)]
                for row in rows]
    kernel = [(-1) ** i * det(centered[:i] + centered[i + 1:])
              for i in range(r + 1)]
    assert any(kernel)
    assert all(sum(kernel[i] * centered[i][j] for i in range(r + 1)) == 0
               for j in range(r))
    assert all(abs(value) <= ell * prod(widths) for value in kernel)
    return True


def prod(values):
    result = 1
    for value in values:
        result *= value
    return result


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def power(z, n):
    out = (1, 0)
    for _ in range(n):
        out = mul(out, z)
    return out


def gaussian_phase_audit(rows, units):
    rows, widths = normalize(rows)
    r = len(widths)
    matrix = [[value - base for value, base in zip(row, rows[0])]
              for row in rows[1:]]
    determinant = det(matrix)
    if not determinant:
        return False
    gaussian_primes = [(2, 1), (3, 2), (4, 1), (5, 2)][:r]
    primes = [a * a + b * b for a, b in gaussian_primes]
    angles = [atan2(b, a) for a, b in gaussian_primes]
    common = (7, 11)
    points = []
    for row, unit in zip(rows, units):
        point = mul(common, power((0, 1), unit))
        for value, width, gaussian in zip(row, widths, gaussian_primes):
            point = mul(point, power(gaussian, value))
            point = mul(point, power((gaussian[0], -gaussian[1]), width - value))
        points.append(point)
    norm = prod(p ** e for p, e in zip(primes, widths))
    norm *= common[0] ** 2 + common[1] ** 2
    assert all(a * a + b * b == norm for a, b in points)
    w = log(norm)
    original_angles = [atan2(b, a) for a, b in points]
    x = [angle - original_angles[0] for angle in original_angles[1:]]
    t = [unit - units[0] for unit in units[1:]]
    lifts = [round((2 * sum(a * b for a, b in zip(row, angles))
                    + pi * unit / 2 - error) / (2 * pi))
             for row, unit, error in zip(matrix, t, x)]
    jmat = adj(matrix)
    for index, row in enumerate(jmat):
        m = 4 * sum(a * b for a, b in zip(row, lifts))
        m -= sum(a * b for a, b in zip(row, t))
        actual = 4 * determinant * angles[index] - pi * m
        expected = 2 * sum(a * b for a, b in zip(row, x))
        assert abs(actual - expected) <= 2e-10 * (1 + abs(expected))
        h, ell = max(widths), factorial(r)
        # delta is an actual containing angular width for these chosen lifts.
        delta = max(original_angles) - min(original_angles)
        assert abs(actual) <= 2 * ell * h ** (r - 1) * delta + 1e-8
        # C=delta*R^(1/2) makes the radius form of the same identity explicit.
        c = delta * norm ** 0.25
        assert abs(actual) <= 2 * ell * c * h ** (r - 1) * norm ** -0.25 + 1e-8
        assert max(2 * abs(determinant), abs(m)) <= (2 + 2 * c) * ell * h ** r
        assert w >= widths[index] * log(primes[index])
    return True


def prime(n):
    return n >= 2 and all(n % divisor for divisor in range(2, isqrt(n) + 1))


def growing_kernel_audit():
    primes = [1009, 1013, 1021]
    assert all(prime(p) and p % 4 == 1 for p in primes)
    count = 0
    for e in [3, 5, 7, 11, 31, 101, 1001]:
        a = (e - 1) // 2
        rows = [[a, a, a], [0, e, e], [e, 0, e], [e, e, 0]]
        assert audit_matrix(rows)
        centered = [[2 * value - e for value in row] for row in rows]
        cofs = [(-1) ** i * det(centered[:i] + centered[i + 1:])
                for i in range(4)]
        divisor = 0
        for value in cofs:
            divisor = gcd(divisor, value)
        if cofs[0] < 0:
            divisor = -divisor
        assert [value // divisor for value in cofs] == [e, 1, 1, 1]
        n, q = [e, 1, 1, 1], e + 3
        assert all(2 * sum(n[i] * rows[i][j] for i in range(4)) == q * e
                   for j in range(3))
        # Exact squared comparison for d_ij > W/2+log(16).
        total_norm = prod(primes) ** e
        for first, second in combinations(rows, 2):
            pair_norm = prod(p ** abs(x - y)
                             for p, x, y in zip(primes, first, second))
            assert pair_norm ** 2 > 256 * total_norm
        # Strictly negative inner products, with denominator e cleared.
        for first, second in combinations(centered, 2):
            assert sum(x * y * log(p) for p, x, y in zip(primes, first, second)) < 0
        count += 1
    return count


def main():
    rng = Random(20260915)
    matrices = 0
    for r in [1, 2]:
        for h in [1, 2, 3]:
            grid = list(product(range(h + 1), repeat=r))
            for rows in combinations(grid, r + 1):
                matrices += audit_matrix([list(row) for row in rows])
    for r in [3, 4]:
        for _ in range(300):
            h = rng.randint(1, 7)
            rows = [[rng.randint(0, h) for _ in range(r)] for _ in range(r + 1)]
            matrices += audit_matrix(rows)
    phases = 0
    for r in [1, 2, 3, 4]:
        for _ in range(80):
            rows = [[rng.randint(0, 3) for _ in range(r)] for _ in range(r + 1)]
            units = [rng.randrange(4) for _ in range(r + 1)]
            phases += gaussian_phase_audit(rows, units)
    kernels = growing_kernel_audit()
    # The coarse numerical specialization used for Matveev's constant.
    assert 2 ** 32 * 4 * 8 * log(2 * 2.71828182846) < 2 ** 40
    print(f"PASS: {matrices} full-rank exact matrices; {phases} actual Gaussian "
          f"phase inversions; {kernels} growing positive kernels.")
    print("Matveev and Roth are cited theorem inputs, not finite-search conclusions.")


if __name__ == "__main__":
    main()
