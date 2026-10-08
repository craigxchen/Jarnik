"""Exact checks for the odd balanced-character quadratic-twist obstruction."""
from itertools import combinations
from math import gcd, isqrt


def rank_f2(rows):
    basis = {}
    for value in rows:
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    return len(basis)


def mul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def conj(z):
    return (z[0], -z[1])


def norm(z):
    return z[0] ** 2 + z[1] ** 2


def power(z, n):
    answer = (1, 0)
    for _ in range(n):
        answer = mul(answer, z)
    return answer


cuts = [c for m in range(1, 5) for c in combinations(range(8), m)
        if m < 4 or 0 in c]
assert [sum(len(c) == m for c in cuts) for m in range(1, 5)] == [8, 28, 56, 35]
assert len(cuts) == 127
signs = [[1 if i in c else -1 for c in cuts] for i in range(8)]
for i in range(8):
    for j in range(8):
        dot = sum(a * b for a, b in zip(signs[i], signs[j]))
        assert dot == 128 * (i == j) - 1
        if i != j:
            assert sum(a != b for a, b in zip(signs[i], signs[j])) == 64
assert 2 * 64 > 127
assert 7 * 8 + 2 * 28 - 56 - 35 == 21

lambdas = [tuple(1 if i in c else -1 for i in range(8))
           for c in combinations(range(8), 4) if 0 in c]
assert len(lambdas) == 35
squareclasses = []
triple_classes = []
for lam in lambdas:
    row = triple_row = 0
    for j, cut in enumerate(cuts):
        exponent = sum(lam[i] for i in cut)
        if len(cut) % 2:
            assert exponent % 2
            bit = 1 << (2 * j + (exponent < 0))
            row |= bit
            if len(cut) == 3:
                triple_row |= bit
        else:
            assert exponent % 2 == 0
    squareclasses.append(row)
    triple_classes.append(triple_row)
assert rank_f2(squareclasses) == rank_f2(triple_classes) == 21
assert rank_f2([row ^ squareclasses[0] for row in squareclasses]) == 20
quadratics = []
for a, b in combinations(range(1, 8), 2):
    quadratics.append(sum(1 << k for k, lam in enumerate(lambdas)
                          if lam[a] == lam[b] == -1))
assert rank_f2(quadratics) == 20
assert rank_f2(quadratics + [(1 << 35) - 1]) == 21
print('PASS: positive 127-cut profile, strict pair separation, full Gram rank, and squareclass rank 21.')


def is_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


blocks = []
p = 5
while len(blocks) < len(cuts):
    if is_prime(p):
        for y in range(1, isqrt(p) + 1):
            x = isqrt(p - y * y)
            if x > y and x * x + y * y == p:
                blocks.append((x, y))
                break
    p += 4
Z = []
for row in signs:
    z = (1, 0)
    for s, block in zip(row, blocks):
        z = mul(z, block if s > 0 else conj(block))
    Z.append(z)
assert len(set(Z)) == 8
assert len({norm(z) for z in Z}) == 1
t = 1
for cut, block in zip(cuts, blocks):
    if len(cut) % 2:
        t *= norm(block)
for lam in lambdas:
    A = (1, 0)
    for cut, block in zip(cuts, blocks):
        exponent = sum(lam[i] for i in cut)
        A = mul(A, power(block if exponent >= 0 else conj(block), abs(exponent)))
    assert gcd(*A) == 1
    assert norm(A) % t == 0
    h = isqrt(norm(A) // t)
    assert norm(A) == t * h * h
    plus = minus = (1, 0)
    for sign, z in zip(lam, Z):
        if sign > 0:
            plus = mul(plus, z)
        else:
            minus = mul(minus, z)
    assert mul(plus, conj(A)) == mul(minus, A)
print('PASS: all 35 literal primitive half-character identities and the common rational norm twist.')

u, y = 1, 0
last_norm = 0
for n in range(1, 36):
    u, y = 9 * u + 20 * y, 4 * u + 9 * y
    x = u - 2 * y
    B = (x, y)
    A = mul((2, 1), mul(B, B))
    assert u * u - 5 * y * y == 1
    assert A[1] == 1 and A[0] < 0
    assert gcd(*A) == 1
    assert norm(A) == 5 * norm(B) ** 2
    assert norm(A) > last_norm
    last_norm = norm(A)
print('PASS: 35 exact fixed-twist Pell examples; the imaginary coordinate is identically one.')
