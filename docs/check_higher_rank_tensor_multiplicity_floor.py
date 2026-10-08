"""Exact checks for the tensor multiplicity floor and width-two sharpness.

The arbitrary-rank tensor-kernel proof is in
higher_rank_tensor_multiplicity_floor.md. The arithmetic fixtures below
are actual primitive Gaussian rows, with one prescribed oriented core.
They do not assert short-arc geometry or a full near-uniform profile.
"""

from itertools import combinations
from math import comb, gcd


def multiplicity_count(n, d):
    q = 2 ** d
    return sum(comb(n, f) * comb(n - f, e) * (q - 2) ** (n - f - e)
               * min(f, e)
               for f in range(n + 1) for e in range(n - f + 1))


def direct_count(n, d):
    answer = 0
    block_mask = (1 << d) - 1
    for subset in range(1 << (n * d)):
        occupancies = [(subset >> (d * a)) & block_mask for a in range(n)]
        answer += min(occupancies.count(block_mask), occupancies.count(0))
    return answer


def matching_pairs(labels):
    if not labels:
        yield ()
        return
    first = labels[0]
    for j in labels[1:]:
        rest = tuple(i for i in labels[1:] if i != j)
        for matching in matching_pairs(rest):
            yield ((first, j),) + matching


def gaussian_multiply(first, second):
    a, b = first
    c, d = second
    return a * c - b * d, a * d + b * c


def valuation(value, prime):
    assert value != 0
    order = 0
    while value % prime == 0:
        value //= prime
        order += 1
    return order


def check_actual_rows(n, depth):
    # pi=10+i has norm 101. These even slopes avoid both primes over 101.
    slopes = [2, 4, 6, 8, 12, 14, 16, 18][:2 * n]
    assert len(slopes) == 2 * n
    core = (1, 0)
    for _ in range(depth):
        core = gaussian_multiply(core, (10, 1))
    matchings = list(matching_pairs(tuple(range(2 * n))))
    for mask in range(1 << (2 * n)):
        T = {i for i in range(2 * n) if mask & (1 << i)}
        rows = [gaussian_multiply(core, (1, slope)) if i in T else (1, slope)
                for i, slope in enumerate(slopes)]
        # Odd norm and gcd of integer coordinates one imply Gaussian primitivity.
        assert all(gcd(x, y) == 1 and (x * x + y * y) % 2 for x, y in rows)
        brackets = {}
        for i, j in combinations(range(2 * n), 2):
            x, y = rows[i]
            xx, yy = rows[j]
            bracket = x * yy - y * xx
            assert bracket != 0
            brackets[i, j] = bracket
            assert valuation(bracket, 101) == (depth if i in T and j in T else 0)
        content = 0
        for matching in matchings:
            product = 1
            for pair in matching:
                product *= brackets[pair]
            content = gcd(content, product)
        assert valuation(content, 101) == depth * max(0, len(T) - n)
        selected = 1
        full = empty = 0
        for a in range(n):
            pair = (2 * a, 2 * a + 1)
            selected *= brackets[pair]
            size = sum(i in T for i in pair)
            full += size == 2
            empty += size == 0
        assert selected % content == 0
        assert valuation(selected // content, 101) == depth * min(full, empty)


def main():
    for n, d in ((2, 2), (2, 5), (3, 3), (4, 2), (4, 3)):
        assert multiplicity_count(n, d) == direct_count(n, d)
    for n in range(2, 51):
        expected = n * 4 ** (n - 1) - n * comb(2 * n, n) // 2
        assert multiplicity_count(n, 2) == expected
    for d in range(2, 21):
        assert multiplicity_count(2, d) == 2
        union = 2 ** (4 * d) - 2 * (2 ** d - 1) ** 4 + (2 ** d - 2) ** 4
        assert multiplicity_count(4, d) == union + 6
    for n in (2, 3, 4):
        for depth in (1, 2):
            check_actual_rows(n, depth)
    print('PASS: tensor multiplicity counts and all-rank width-two formula.')
    print('PASS: exact primitive coefficient orders for 672 Gaussian core fixtures.')


if __name__ == '__main__':
    main()
