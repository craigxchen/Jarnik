"""Exact rank/degree and first-syzygy checks for the EN kernel resolution.

The general exactness and codimension-two extension are prose arguments.
"""

from functools import lru_cache
from math import comb, factorial
from random import Random

from check_higher_rank_terminal_fusion_grid import case


def tableau_hook(shape):
    shape = tuple(x for x in shape if x)
    hooks = 1
    for i, width in enumerate(shape):
        for j in range(width):
            hooks *= width - j + sum(other > j for other in shape[i + 1:])
    return factorial(sum(shape)) // hooks


@lru_cache(None)
def tableau_branch(shape):
    if not shape:
        return 1
    answer = 0
    for i in range(len(shape)):
        if i + 1 == len(shape) or shape[i] > shape[i + 1]:
            reduced = list(shape)
            reduced[i] -= 1
            answer += tableau_branch(tuple(x for x in reduced if x))
    return answer


def resolution(n, d):
    e = d - 2
    m = n * d
    coefficients = [0] * (2 * e + 1)
    for a in range(e + 1):
        for b in range(a + 1):
            shape = (e,) * (n - 2) + (e - b, e - a)
            dimension = tableau_hook(shape)
            if n <= 5 and d <= 6:
                assert dimension == tableau_branch(tuple(x for x in shape if x))
            coefficients[a + b] += (a - b + 1) * dimension
    ranks = [comb(m, 2 * n + j) * value
             for j, value in enumerate(coefficients)]
    degrees = [comb(m - 1, 2 * n + j - 1) * value
               for j, value in enumerate(coefficients)]
    h = sum((-1) ** j * value for j, value in enumerate(ranks))
    delta = sum((-1) ** j * value for j, value in enumerate(degrees))
    return coefficients, ranks, h, delta


def det_mod(matrix, prime):
    matrix = [[value % prime for value in row] for row in matrix]
    result = 1
    for column in range(len(matrix)):
        pivot = next((i for i in range(column, len(matrix))
                      if matrix[i][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            matrix[pivot], matrix[column] = matrix[column], matrix[pivot]
            result = -result
        value = matrix[column][column]
        result = result * value % prime
        inverse = pow(value, -1, prime)
        for i in range(column + 1, len(matrix)):
            multiplier = matrix[i][column] * inverse % prime
            for j in range(column + 1, len(matrix)):
                matrix[i][j] = (matrix[i][j]
                                 - multiplier * matrix[column][j]) % prime
    return result % prime


def laplace_fixtures():
    prime = 1009
    random = Random(20261001)
    tested = 0
    for n in range(2, 5):
        for fixture in range(5):
            matrix = []
            for i in range(2 * n + 1):
                p = 1 + (i + 1) * (fixture + 1)
                y = (i + 2) ** 2 + fixture
                v = [random.randrange(1, prime) for _ in range(n)]
                matrix.append([p * x for x in v] + [y * x for x in v])
            cofactors = [(-1) ** i * det_mod(matrix[:i] + matrix[i + 1:], prime)
                         for i in range(2 * n + 1)]
            assert any(cofactors), 'A zero cofactor vector would make this vacuous'
            for column in range(2 * n):
                assert sum(row[column] * cofactor
                           for row, cofactor in zip(matrix, cofactors)) % prime == 0
                tested += 1
    return tested


def main():
    for n in range(2, 9):
        for d in range(2, 13):
            _, _, h, delta = resolution(n, d)
            fusion = case(n, d)
            assert (h, delta) == (fusion['h'], fusion['delta'])
    coefficients, ranks, h, delta = resolution(3, 4)
    assert coefficients == [5, 10, 9, 4, 1]
    assert ranks == [4620, 7920, 4455, 880, 66]
    assert (h, delta) == (341, 55)
    print('PASS: 77 EN rank/degree formulas agree with independent fusion data.')
    print('PASS: hook dimensions agree with independent tableau branching.')
    print('Twelve-row coefficients:', coefficients, 'ranks:', ranks)
    print(f'PASS: {laplace_fixtures()} exact first-syzygy identities.')


if __name__ == '__main__':
    main()
