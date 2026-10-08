"""Exhaustive four-term balanced relations and their sine-factor heights."""

from collections import Counter, defaultdict
from itertools import combinations, product
from fractions import Fraction

ROWS = range(8)
L = [v for v in product((-1, 1), repeat=8) if sum(v) == 0 and v[0] == 1]
assert len(L) == 35

def add(a, b):
    return tuple(x + y for x, y in zip(a, b))

def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))

def scale2(a):
    assert all(x % 2 == 0 for x in a)
    return tuple(x // 2 for x in a)

def dot(a, b):
    return sum(x * y for x, y in zip(a, b))

def support(a):
    return sum(x != 0 for x in a)

pair_sum = defaultdict(list)
for i, lam in enumerate(L):
    for j in range(i, len(L)):
        pair_sum[add(lam, L[j])].append((i, j))
assert Counter(map(len, pair_sum.values())) == Counter({1: 245, 3: 105, 10: 7})

patterns = []
for pairs in pair_sum.values():
    distinct = [pair for pair in pairs if pair[0] != pair[1]]
    for left, right in combinations(distinct, 2):
        i, j = left
        k, m = right
        lam, mu, nu, rho = L[i], L[j], L[k], L[m]
        assert add(lam, mu) == add(nu, rho)
        p = scale2(sub(lam, mu))
        q = scale2(sub(nu, rho))
        r = scale2(add(p, q))
        s = scale2(sub(p, q))
        assert all(sum(v) == 0 for v in (p, q, r, s))
        assert all(v in (-1, 0, 1) for v in r + s)
        assert len(set((i, j, k, m))) == 4
        patterns.append((p, q, r, s))
assert len(patterns) == 630

classification = Counter((support(p), support(q), dot(p, q), support(r), support(s))
                         for p, q, r, s in patterns)
assert classification == Counter({
    (4, 4, 0, 2, 2): 315,
    (6, 6, 2, 4, 2): 210,
    (6, 6, -2, 2, 4): 105,
})
print("PASS: 630 relations and exact (315,210,105) support classification.")

columns = [c for c in product((-1, 1), repeat=8) if abs(sum(c)) < 8]
assert len(columns) == 254

factors = {f for _, _, r, s in patterns for f in (r, s)}
for representative in factors:
    size = support(representative)
    by_cut_size = {}
    for t in range(1, 5):
        hist = Counter(Fraction(abs(dot(representative, c)), 2)
                       for c in columns if min(c.count(1), c.count(-1)) == t)
        by_cut_size[t] = tuple(hist[k] for k in (0, 1, 2))
    expected = {
        2: {1: (12, 4, 0), 2: (32, 24, 0), 3: (52, 60, 0), 4: (30, 40, 0)},
        4: {1: (8, 8, 0), 2: (20, 32, 4), 3: (40, 56, 16), 4: (26, 32, 12)},
    }[size]
    assert by_cut_size == expected
    full = Counter(Fraction(abs(dot(representative, c)), 2) for c in columns)
    expected_full = {0: (126, 94)[size == 4],
                     1: (128, 128)[size == 4],
                     2: (0, 32)[size == 4]}
    assert full == Counter({k: v for k, v in expected_full.items() if v})
    assert sum(Fraction(abs(dot(representative, c)), 2) for c in columns) == (128, 192)[size == 4]
print(f"PASS: all {len(factors)} occurring factors and 254 columns have exact height histograms.")

# The sine factorization is the cosine-difference identity in vector form.
# Check its four exponent vectors exactly, for every relation.
for p, q, r, s in patterns:
    assert add(r, s) == p and sub(r, s) == q
print("PASS: every relation factors through the two integral sine vectors.")

# Formal equal-weight 127-cut profile: divide the oriented sign-column sum by
# two to pass to one representative of each complementary cut.
# Use a direct canonical representative for each complementary cut.
unoriented = [c for c in columns if c[0] == 1]
assert len(unoriented) == 127
for size, target in ((2, 64), (4, 96)):
    f = next(f for _, _, r, s in patterns for f in (r, s) if support(f) == size)
    assert sum(Fraction(abs(dot(f, c)), 2) for c in unoriented) == target
print("PASS: equal-weight leading heights are 64h and 96h; additive totals are 128h or 160h.")
