"""Exact Gaussian content and scalar-height checks for torus relations."""
from itertools import product
from random import Random


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def mul(z, w):
    return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])


def conj(z):
    return (z[0], -z[1])


def power(z, n):
    out = (1, 0)
    while n:
        if n & 1:
            out = mul(out, z)
        z = mul(z, z)
        n //= 2
    return out


def oriented_power(z, n):
    return power(z if n >= 0 else conj(z), abs(n))


characters = [v for v in product((-1, 1), repeat=8)
              if sum(v) == 0 and v[0] == 1]
blocks = [(2, 1), (3, 2), (4, 1), (5, 2)]
primes = [5, 13, 17, 29]
rng = Random(20260921)
nontrivial_content = 0
for trial in range(120):
    allocations = [[rng.randrange(8) for _ in range(8)] for _ in blocks]
    cs = [[dot(lam, a) for a in allocations] for lam in characters]
    A = []
    for c in cs:
        z = (1, 0)
        for h, exponent in zip(blocks, c):
            z = mul(z, oriented_power(h, exponent))
        A.append(z)
    n = [0]*35
    for index in rng.sample(range(35), 8):
        n[index] = rng.choice((-2, -1, 1, 2))
    row = [sum(n[k]*characters[k][i] for k in range(35)) for i in range(8)]
    assert sum(row) == 0
    G = (1, 0)
    for z, exponent in zip(A, n):
        G = mul(G, oriented_power(z, exponent))
    C = 1
    primitive = (1, 0)
    for j, (h, p, a) in enumerate(zip(blocks, primes, allocations)):
        q = dot(row, a)
        assert q == sum(n[k]*cs[k][j] for k in range(35))
        total = sum(abs(n[k])*abs(cs[k][j]) for k in range(35))
        assert total >= abs(q) and (total-abs(q)) % 2 == 0
        C *= p**((total-abs(q))//2)
        primitive = mul(primitive, oriented_power(h, q))
    assert G == (C*primitive[0], C*primitive[1])
    nontrivial_content += C > 1
assert nontrivial_content > 0

allocation = (0, 0, 5, 4, 1, 0, 0, 0)
negative_sets = ({0, 1, 2, 3}, {0, 1, 4, 5},
                 {0, 1, 2, 4}, {0, 1, 3, 5})
lam = [tuple(-1 if i in s else 1 for i in range(8)) for s in negative_sets]
assert all(lam[0][i]+lam[1][i] == lam[2][i]+lam[3][i] for i in range(8))
exponents = [dot(v, allocation) for v in lam]
assert exponents == [-8, 8, -2, 2]
values = [oriented_power((2, 1), e) for e in exponents]
assert mul(values[0], values[1]) == (5**8, 0)
assert mul(values[2], values[3]) == (5**2, 0)
print('PASS: 120 exact nested Gaussian content identities and the p^6 exchange ratio.')

count = 0
minimum = None
for first in product(range(-2, 3), repeat=7):
    if not any(first):
        continue
    row = first + (-sum(first),)
    subset_sums = [0]
    for coefficient in row:
        subset_sums += [x+coefficient for x in subset_sums]
    total = sum(abs(x) for x in subset_sums)
    assert total % 2 == 0
    cut_height = total//2
    assert cut_height >= 64*max(abs(x) for x in row)
    minimum = cut_height if minimum is None else min(minimum, cut_height)
    count += 1
assert count == 78124 and minimum == 64
print('PASS: 78,124 bounded integer row vectors; minimum full-cut scalar height 64.')
