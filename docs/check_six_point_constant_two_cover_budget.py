"""Verify finite combinatorial bounds for the constant-two-cover audit.

python3 docs/check_six_point_constant_two_cover_budget.py
The modular polynomial certificates are checked separately by
check_six_point_branch_two_torsion.py.
"""
from math import factorial, prod, comb
from fractions import Fraction


def partitions(n, minimum=1):
    if n == 0:
        yield ()
    for first in range(minimum, n + 1):
        for tail in partitions(n - first, first):
            yield (first,) + tail


def main():
    for m in range(2, 11):
        products = []
        for counts in partitions(12):
            k = len(counts)
            if k >= m + 2 or k == m + 1 and all(c % 2 == 0 for c in counts):
                products.append(prod(factorial(c) for c in counts))
        assert max(products) <= factorial(11 - m)
        if m >= 3:
            degree_lower = Fraction(factorial(12), 8 * factorial(11 - m))
            assert degree_lower / 2 > 40 * 2**m
            print(f'PASS: rank={m}, archimedean lower={degree_lower / 2}, '
                  f'geometric pole upper={40 * 2**m}')
    maximum_stabilizer = 0
    for a in range(13):
        for b in range(13 - a):
            for c in range(13 - a - b):
                counts = (a, b, c, 12 - a - b - c)
                if sum(n > 0 for n in counts) < 3:
                    continue
                if (counts[1] + counts[3]) % 2 or (counts[2] + counts[3]) % 2:
                    continue
                translations = [t for t in range(4)
                                if all(counts[x] == counts[x ^ t] for x in range(4))]
                stabilizer = len(translations) * prod(factorial(n) for n in counts)
                maximum_stabilizer = max(maximum_stabilizer, stabilizer)
    assert maximum_stabilizer == factorial(9)
    degree_lower = factorial(12) // (2 * maximum_stabilizer)
    assert degree_lower == 660 and degree_lower / 2 > 40 * 4
    print('PASS: rank=2 exact maximum S12 pointwise stabilizer=9!; '
          'degree>=660, archimedean>=330>160.')
    orbit_sizes = (comb(12, 2), comb(12, 4), comb(12, 6) // 2)
    assert orbit_sizes == (66, 495, 462)
    assert min(orbit_sizes) / 2 > 20
    print('PASS: rank=1 orbit sizes=(66,495,462); '
          'archimedean>=33e>20e pole orbits after extension degree e.')


if __name__ == '__main__':
    main()
