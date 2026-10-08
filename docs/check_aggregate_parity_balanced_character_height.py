"""Exact aggregate-parity heights and the anchored ratio-averaging correction."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import gcd
from random import Random

L = [v for v in product((-1, 1), repeat=8) if sum(v) == 0 and v[0] == 1]
assert len(L) == 35

def dot(a, b):
    return sum(x*y for x, y in zip(a, b))

expected = {
    1: Counter({Fraction(1, 2): 35}),
    2: Counter({0: 20, 1: 15}),
    3: Counter({Fraction(1, 2): 30, Fraction(3, 2): 5}),
    4: Counter({0: 18, 1: 16, 2: 1}),
}
columns = [c for c in product((-1, 1), repeat=8) if abs(sum(c)) < 8]
for c in columns:
    r = min(c.count(1), c.count(-1))
    values = [Fraction(abs(dot(lam, c)), 4) for lam in L]
    assert Counter(values) == expected[r]
assert len(columns) == 254
assert [sum(v*n for v, n in expected[r].items()) for r in range(1, 5)] == [Fraction(35, 2), 15, Fraction(45, 2), 18]

anchored = []
folded = []
pairs = Counter()
for lam, mu in combinations(L, 2):
    distance = sum(x != y for x, y in zip(lam, mu))
    if distance not in (2, 6):
        continue
    v = tuple(x-y if distance == 2 else x+y for x, y in zip(lam, mu))
    support = tuple(i for i, x in enumerate(v) if x)
    assert len(support) == 2 and sorted(v[i] for i in support) == [-2, 2]
    pairs[support] += 1
    folded.append(v)
    if distance == 2:
        anchored.append(v)
assert len(anchored) == 210 and len(folded) == 280
assert pairs == Counter({pair: 10 for pair in combinations(range(8), 2)})
for c in columns:
    k = sum(x == -1 for x in c[1:])
    anchored_total = sum(Fraction(abs(dot(v, c)), 4) for v in anchored)
    assert anchored_total == 10*k*(7-k)
    r = min(c.count(1), c.count(-1))
    assert sum(Fraction(abs(dot(v, c)), 4) for v in folded) == 10*r*(8-r)
print('PASS: all 254 cut tables; 210 anchored edges and 280 folded edges = ten copies of all 28 row pairs.')

# Exact source-prime cancellation, including odd threshold layers.
for x in range(8):
    for y in range(8):
        a = (0, x, x, x, x, x, x, x+y)
        actual = sum(Fraction(abs(dot(lam, a)), 2) for lam in L)
        upper = Fraction(35, 2)*(x+y)
        assert upper-actual == 15*min(x, y)
        assert all(dot(lam, a) % 2 == (x+y) % 2 for lam in L)
a = (0, 0, 0, 0, 0, 0, 0, 2)
assert sum(Fraction(abs(dot(lam, a)), 2) for lam in L) == 35

# Mixed nested allocations: the total allocation is even at each prime,
# while individual odd threshold layers occur. Repeated blocks cancel.
def mul(z, w):
    return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])

def conj(z):
    return (z[0], -z[1])

def norm(z):
    return z[0]*z[0]+z[1]*z[1]

def power(z, n):
    out = (1, 0)
    for _ in range(n):
        out = mul(out, z)
    return out

blocks = [(2, 1), (3, 2), (4, 1), (5, 2)]
assert [norm(z) for z in blocks] == [5, 13, 17, 29]
rng = Random(20260921)
units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
cases = cancellation_cases = odd_layers = 0
for trial in range(24):
    allocations = []
    for _ in blocks:
        a = [rng.randrange(7) for _ in range(8)]
        if sum(a) % 2:
            a[0] += 1
        assert sum(a) % 2 == 0
        allocations.append(a)
        odd_layers += sum(sum(t >= k for t in a) % 2 for k in range(1, max(a)+1))
    exps = [rng.randrange(4) for _ in range(8)]
    Z = []
    for i in range(8):
        z = mul((2, 1), units[exps[i]])
        for h, a in zip(blocks, allocations):
            z = mul(z, mul(power(h, a[i]), power(conj(h), max(a)-a[i])))
        Z.append(z)
    assert len({norm(z) for z in Z}) == 1
    for lam in L:
        beta = (1, 0)
        for h, a in zip(blocks, allocations):
            n = dot(lam, a)
            assert n % 2 == 0
            e = n//2
            upper = sum(Fraction(abs(dot(lam, [1 if x >= k else -1 for x in a])), 4)
                        for k in range(1, max(a)+1))
            assert upper >= abs(e)
            cancellation_cases += upper > abs(e)
            beta = mul(beta, power(h if e >= 0 else conj(h), abs(e)))
        assert gcd(abs(beta[0]), abs(beta[1])) == 1 and norm(beta) % 2 == 1
        plus = minus = (1, 0)
        for l, z in zip(lam, Z):
            if l > 0:
                plus = mul(plus, z)
            else:
                minus = mul(minus, z)
        u = units[dot(lam, exps) % 4]
        assert mul(plus, power(conj(beta), 2)) == mul(mul(u, minus), power(beta, 2))
        cases += 1
assert cases == 840 and odd_layers > 0 and cancellation_cases > 0
print(f'PASS: singleton cancellation identity; {cases} mixed nested Gaussian identities with {odd_layers} odd layers.')
