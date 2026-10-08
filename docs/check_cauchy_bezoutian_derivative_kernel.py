"""Exact Cauchy/Bezoutian factorization and derivative-kernel audit."""

from fractions import Fraction as F
from itertools import combinations
from math import gcd, prod

from check_integer_cotangent_positive_bezoutian import (
    bezoutian, determinant, evaluate,
)


def multiply(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)]
            for row in a]


def rank(matrix):
    a = [[F(x) for x in row] for row in matrix]
    pivot = 0
    for col in range(len(a[0])):
        index = next((i for i in range(pivot, len(a)) if a[i][col]), None)
        if index is None:
            continue
        a[pivot], a[index] = a[index], a[pivot]
        scale = a[pivot][col]
        a[pivot] = [v/scale for v in a[pivot]]
        for i in range(len(a)):
            if i != pivot:
                scale = a[i][col]
                a[i] = [x-scale*y for x, y in zip(a[i], a[pivot])]
        pivot += 1
        if pivot == len(a):
            break
    return pivot


def rational_determinant(matrix):
    a = [[F(x) for x in row] for row in matrix]
    result = F(1)
    for col in range(len(a)):
        index = next((i for i in range(col, len(a)) if a[i][col]), None)
        if index is None:
            return F(0)
        if index != col:
            a[index], a[col] = a[col], a[index]
            result = -result
        result *= a[col][col]
        for i in range(col+1, len(a)):
            scale = a[i][col]/a[col][col]
            a[i] = [x-scale*y for x, y in zip(a[i], a[col])]
    return result


def check(xs):
    k, r = len(xs), len(xs)//2
    assert len(set(xs)) == k and all(x > 0 for x in xs)
    _, _, o, b = bezoutian(xs)
    os = [prod(x+y for j, y in enumerate(xs) if i != j)
          for i, x in enumerate(xs)]
    ds = [prod(y-x for j, y in enumerate(xs) if i != j)
          for i, x in enumerate(xs)]
    assert os == [evaluate(o, x*x) for x in xs]
    ws = [F(d, 2*x*oi) for x, oi, d in zip(xs, os, ds)]
    c = [[F(1, x+y) for y in xs] for x in xs]
    v = [[F(x**(2*a), oi) for a in range(r)] for x, oi in zip(xs, os)]
    gram = multiply(multiply(v, b), list(zip(*v)))
    minus = [[c[i][j]-(ws[i] if i == j else 0) for j in range(k)]
             for i in range(k)]
    plus = [[c[i][j]+(ws[i] if i == j else 0) for j in range(k)]
            for i in range(k)]
    assert gram == minus
    assert rank(minus) == r and rank(plus) == k-r
    kernel = [[F(x**(2*q), (-1)**(k-1)*d) for x, d in zip(xs, ds)]
              for q in range(k-r)]
    assert rank(kernel) == k-r
    assert not any(any(row) for row in multiply(minus, list(zip(*kernel))))
    rho = [(-1)**(k-1)/w for w in ws]
    inverse = [[rho[i]*c[i][j]*rho[j] for j in range(k)] for i in range(k)]
    assert multiply(c, inverse) == [[int(i == j) for j in range(k)] for i in range(k)]
    assert [[ws[i]*inverse[i][j]*ws[j] for j in range(k)] for i in range(k)] == c
    positives = [i for i, w in enumerate(ws) if w > 0]
    for size in range(1, len(positives)+1):
        for subset in combinations(positives, size):
            lhs = rational_determinant([[c[i][j] for j in subset] for i in subset])
            lhs /= prod(ws[i] for i in subset)
            rhs = prod(F(xs[i]+xs[j], abs(xs[i]-xs[j]))
                       for i in subset for j in range(k) if j not in subset)
            assert lhs == rhs >= 1
    return minus


def main():
    cases = 0
    for k in range(2, 13):
        for xs in (list(range(1, k+1)), [2**j+1 for j in range(k)],
                   [10000+3*j+j*j for j in range(k)]):
            check(xs)
            cases += 1
    xs = [3456, 3723, 7184]
    matrix = check(xs)
    assert bezoutian(xs)[3] == [[833126395920]]
    assert multiply(matrix, [[10640], [-10907], [0]]) == [[0], [0], [0]]
    xs = [1146, 1257, 1842, 4104]
    check(xs)
    b = bezoutian(xs)[3]
    assert b == [[531601657320495744, 26730985716], [26730985716, 8349]]
    assert gcd(*(x for row in b for x in row)) == 3
    det = determinant(b)
    assert det == 3723796639619822934000 and det % 125 == 0 and det % 625
    assert all(det % p for p in (13, 17, 29, 37, 73))
    assert b[1][1]//3 == 2783
    print(f'PASS: {cases+2} exact Cauchy/Bezoutian cases; factorization, derivative')
    print('kernels, both complementary ranks, inverse, principal-minor identity, fixtures.')


if __name__ == '__main__':
    main()
