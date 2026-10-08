"""Exact checks for positive_hankel_full_profile_budget.md.

Standard-library only. Finite checks supplement the general proofs.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt, lcm, prod


def popcount(value):
    return bin(value).count("1")


def gaussian_mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def gaussian_norm(z):
    return z[0] * z[0] + z[1] * z[1]


def gaussian_prime(p):
    for a in range(1, isqrt(p) + 1):
        b2 = p - a * a
        b = isqrt(b2)
        if b > 0 and b * b == b2:
            return a, b
    raise AssertionError("not a split prime: %s" % p)


def determinant(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    n = len(a)
    answer = Fraction(1)
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            answer = -answer
        value = a[j][j]
        answer *= value
        for i in range(j + 1, n):
            q = a[i][j] / value
            for k in range(j + 1, n):
                a[i][k] -= q * a[j][k]
            a[i][j] = 0
    return answer


def selection_masks(m, r):
    return [sum(1 << i for i in s) for s in combinations(range(m), r)]


def block_exponent(t_mask, s_mask):
    s = popcount(t_mask)
    t = popcount(t_mask & s_mask)
    return s - t + t * (t - 1)


def check_combinatorial_budget():
    cases = 0
    for m in range(2, 10):
        for r in range(1, m + 1):
            selections = selection_masks(m, r)
            blocks = range(1, 1 << m)
            size_exponent = ((m - r) * (1 << (m - 1))
                             + r * (r - 1) * (1 << (m - 2)))
            forced_exponent = sum(
                min(block_exponent(t, s) for s in selections)
                for t in blocks
            )
            assert all(
                sum(block_exponent(t, s) for t in blocks) == size_exponent
                for s in selections
            )
            assert (forced_exponent < size_exponent) == (r < m)
            if m == 4:
                assert forced_exponent == (17, 18, 25, 48)[r - 1]
                assert size_exponent == (24, 24, 32, 48)[r - 1]
            cases += 1
    return cases


def check_full_fifteen_block_fixture():
    primes = (5, 13, 17, 29, 37, 41, 53, 61, 73,
              89, 97, 101, 109, 113, 137)
    factors = {mask: gaussian_prime(p)
               for mask, p in enumerate(primes, 1)}
    assert all(gaussian_norm(factors[mask]) == primes[mask - 1]
               for mask in factors)

    rows = []
    for i in range(4):
        row = (1, 0)
        for mask in range(1, 16):
            if mask & (1 << i):
                row = gaussian_mul(row, factors[mask])
        # A correction overlapping the incident singleton core block.
        row = gaussian_mul(row, factors[1 << i])
        rows.append(row)
    assert all(y and gcd(x, y) == 1 and (x - y) % 2
               for x, y in rows)

    residues = []
    for i, j in combinations(range(4), 2):
        common = prod(primes[mask - 1] for mask in range(1, 16)
                      if mask & (1 << i) and mask & (1 << j))
        delta = rows[i][0] * rows[j][1] - rows[j][0] * rows[i][1]
        assert delta and delta % common == 0
        residues.append(abs(delta // common))
    scale = lcm(*(abs(y) for _, y in rows), *residues)
    cotangents = [scale * x // y for x, y in rows]
    norms = [x * x + scale * scale for x in cotangents]
    assert len(set(cotangents)) == 4

    for mask, p in enumerate(primes, 1):
        incident = [i for i in range(4) if mask & (1 << i)]
        assert all(norms[i] % p == 0 for i in incident)
        assert all((cotangents[i] - cotangents[j]) % p == 0
                   for i, j in combinations(incident, 2))

    checks = 0
    totals = {}
    for r in range(1, 5):
        total = 0
        for s in combinations(range(4), r):
            mask = sum(1 << i for i in s)
            total += (
                prod((cotangents[i] - cotangents[j]) ** 2
                     for i, j in combinations(s, 2))
                * prod(norms[i] for i in range(4) if not mask & (1 << i))
            )
        moment = [
            [sum(Fraction(cotangents[i] ** (a + b), norms[i])
                 for i in range(4)) for b in range(r)]
            for a in range(r)
        ]
        gram_value = prod(norms) * determinant(moment)
        assert gram_value.denominator == 1
        assert gram_value.numerator == total > 0

        selections = selection_masks(4, r)
        guaranteed = prod(
            p ** min(block_exponent(mask, s) for s in selections)
            for mask, p in enumerate(primes, 1)
        )
        assert total % guaranteed == 0
        totals[r] = total
        checks += 1

    # Formula (10): the first surviving residue at singleton core primes.
    for i, expected_d2 in ((1, 12), (2, 6), (3, 1)):
        p = primes[(1 << i) - 1]
        assert scale % p and all(norms[j] % p for j in range(4) if j != i)
        alpha = cotangents[i] % p
        assert (alpha * alpha + scale * scale) % p == 0
        others = tuple(j for j in range(4) if j != i)
        rho = {j: ((cotangents[j] - alpha)
                   * pow(cotangents[j] + alpha, p - 2, p)) % p
               for j in others}
        outside_norms = prod(norms[j] for j in others)
        for r in range(1, 5):
            residue = sum(
                prod((cotangents[j] - cotangents[k]) ** 2
                     for j, k in combinations(subset, 2))
                * prod(rho[j] for j in subset)
                for subset in combinations(others, r - 1)
            )
            assert totals[r] % p == outside_norms * residue % p
        assert totals[2] % p == expected_d2
    return checks


def main():
    budgets = check_combinatorial_budget()
    minors = check_full_fifteen_block_fixture()
    print("PASS: %d exact combinatorial budgets, m=2..9" % budgets)
    print("PASS: %d positive Hankel minors on a primitive full-15-block fixture" % minors)
    print("PASS: all core divisibilities; corrections overlap singleton core factors")
    print("PASS: 12 singleton-prime leading-residue identities, including nonzero D2")


if __name__ == "__main__":
    main()
