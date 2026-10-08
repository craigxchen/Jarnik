"""Exact checks for additive normalized balanced-character content.

The all-order support inequality is proved in the cited multiplicity note;
these tests check its arithmetic normalization and a sharp support example.
"""
from fractions import Fraction
from itertools import combinations, product
from math import gcd, prod
from random import Random


def mul(z, w):
    return z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0]


def conj(z):
    return z[0], -z[1]


def power(z, n):
    out = (1, 0)
    while n:
        if n & 1:
            out = mul(out, z)
        z = mul(z, z)
        n //= 2
    return out


def oriented(z, n):
    return power(z if n >= 0 else conj(z), abs(n))


def dot(v, a):
    return sum(x*y for x, y in zip(v, a))


def gaussian_product(values):
    out = (1, 0)
    for value in values:
        out = mul(out, value)
    return out


characters = [v for v in product((-1, 1), repeat=8) if sum(v) == 0]
blocks = [(2, 1), (3, 2), (4, 1), (5, 2)]
primes = [5, 13, 17, 29]
rng = Random(202609211)
bilinear_count = 0
clearing_count = 0
nontrivial_content = False
odd_twist_seen = False
for trial in range(24):
    allocations = [[rng.randrange(7) for _ in range(8)] for _ in blocks]
    cs = [[dot(lam, a) for a in allocations] for lam in characters]
    values = [gaussian_product(oriented(z, c) for z, c in zip(blocks, row))
              for row in cs]
    parity = [sum(a) % 2 for a in allocations]
    t = prod(p**e for p, e in zip(primes, parity))
    hs = [prod(p**((abs(c)-e)//2)
               for p, c, e in zip(primes, row, parity)) for row in cs]
    odd_twist_seen |= t > 1
    for value, h in zip(values, hs):
        assert value[0]**2+value[1]**2 == t*h*h
    for _ in range(12):
        j, k = rng.sample(range(70), 2)
        sigma = rng.choice((-1, 1))
        qs = [(a+sigma*b)//2 for a, b in zip(cs[j], cs[k])]
        B = gaussian_product(oriented(z, q) for z, q in zip(blocks, qs))
        C = prod(p**((abs(a)+abs(b)-2*abs(q))//2)
                 for p, a, b, q in zip(primes, cs[j], cs[k], qs))
        actual = mul(values[j], values[k] if sigma == 1 else conj(values[k]))
        square = mul(B, B)
        assert actual == (C*square[0], C*square[1])
        x, y = B
        assert actual[0] == C*(x*x-y*y)
        assert actual[1] == 2*C*x*y
        assert C*(x*x+y*y)-actual[0] == 2*C*y*y
        assert gcd(abs(actual[0]), abs(actual[1])) == C
        nontrivial_content |= C > 1
        bilinear_count += 1
    for d in range(1, 5):
        terms = [tuple(rng.randrange(70) for _ in range(d)) for _ in range(6)]
        signed_exponents = [[sum(cs[j][p] for j in term) for p in range(4)]
                            for term in terms]
        # Include conjugates: the supported exponent set is reciprocal.
        m = [max(abs(c[p]) for c in signed_exponents) for p in range(4)]
        kappa = [d*e % 2 for e in parity]
        assert all((bound-e) % 2 == 0 for bound, e in zip(m, kappa))
        Q = prod(p**((bound-e)//2) for p, bound, e in zip(primes, m, kappa))
        td = prod(p**e for p, e in zip(primes, kappa))
        assert td == (t if d % 2 else 1)
        assert Q*Q*td == prod(p**bound for p, bound in zip(primes, m))
        for term, cs_term in zip(terms, signed_exponents):
            numerator = gaussian_product(values[j] for j in term)
            denominator = t**(d//2)*prod(hs[j] for j in term)
            scaled = tuple(Fraction(Q*x, denominator) for x in numerator)
            assert all(x.denominator == 1 for x in scaled)
            direct = gaussian_product(
                mul(power(z, (bound+c)//2), power(conj(z), (bound-c)//2))
                for z, bound, c in zip(blocks, m, cs_term))
            assert scaled == direct
            clearing_count += 1

assert nontrivial_content and odd_twist_seen
print(f'PASS: {bilinear_count} nested-prime bilinears, exact rational contents.')
print(f'PASS: {clearing_count} twist-cleared monomials in degrees 1 through 4.')

cuts = []
for r in range(1, 5):
    cuts += [frozenset(S) for S in combinations(range(8), r)
             if r < 4 or 0 in S]
assert len(cuts) == 127
assert [sum(len(S) == r for S in cuts) for r in range(1, 5)] == [8, 28, 56, 35]
assert sum(max(abs(sum(lam[i] for i in S)) for lam in characters)
           for S in cuts) == 372

# Expansion of product_j(z_2j-z_(2j+1)), a reciprocal balanced support.
matching_support = {}
for choice in product((0, 1), repeat=4):
    selected = frozenset(2*j+choice[j] for j in range(4))
    lam = tuple(1 if i in selected else -1 for i in range(8))
    matching_support[lam] = (-1)**sum(choice)
assert len(matching_support) == 16
assert all(matching_support[tuple(-x for x in lam)] == coefficient
           for lam, coefficient in matching_support.items())
total_twice_width = 0
for S in cuts:
    twice_width = max(abs(sum(lam[i] for i in S)) for lam in matching_support)
    separated = sum((2*j in S) != (2*j+1 in S) for j in range(4))
    assert twice_width == separated
    total_twice_width += twice_width
assert total_twice_width == 256

# Exact Taylor expansion about all ones: coefficients vanish below degree4.
taylor = {}
for lam, coefficient in matching_support.items():
    selected = tuple(i for i, x in enumerate(lam) if x == 1)
    for mask in range(16):
        monomial = tuple(selected[j] for j in range(4) if mask >> j & 1)
        taylor[monomial] = taylor.get(monomial, 0)+coefficient
nonzero_taylor = {monomial: c for monomial, c in taylor.items() if c}
assert len(nonzero_taylor) == 16
assert all(len(monomial) == 4 for monomial in nonzero_taylor)
assert Fraction(total_twice_width, 2) == 32*4
assert Fraction(total_twice_width, 2)-Fraction(127*4, 4) == 1
print('PASS: all 127 cuts; sharp order-four matching support has width128 and margin1.')
