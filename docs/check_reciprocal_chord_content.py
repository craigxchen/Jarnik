"""Exact arithmetic checks for reciprocal_chord_content.md.

These are actual equal-norm Gaussian tuples, not endpoint examples.
All theorem proofs in the note are symbolic and allow unbounded exponents.
"""

from itertools import combinations, permutations, product
from math import atan2, comb, gcd, prod
from random import Random


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def mul(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def conj(a):
    return a[0], -a[1]


def norm(a):
    return a[0] ** 2 + a[1] ** 2


def gp(a, n):
    answer = (1, 0)
    for _ in range(n):
        answer = mul(answer, a)
    return answer


def gprod(values):
    answer = (1, 0)
    for value in values:
        answer = mul(answer, value)
    return answer


def divide(a, b):
    numerator = mul(a, conj(b))
    denominator = norm(b)
    assert all(t % denominator == 0 for t in numerator)
    return tuple(t // denominator for t in numerator)


def ggcd(a, b):
    while b != (0, 0):
        numerator = mul(a, conj(b))
        denominator = norm(b)
        quotient = tuple((2 * t + denominator) // (2 * denominator)
                         for t in numerator)
        remainder = add(a, tuple(-t for t in mul(quotient, b)))
        assert norm(remainder) < norm(b)
        a, b = b, remainder
    return a


def glcm(values):
    answer = (1, 0)
    for value in values:
        answer = mul(divide(answer, ggcd(answer, value)), value)
    return answer


PERMS = list(permutations(range(4)))
SIGNS = [(-1) ** sum(s[i] > s[j] for i in range(4) for j in range(i + 1, 4))
         for s in PERMS]
PRIMES = [(2, 1), (10, 1), (14, 1)]
PNORMS = [norm(p) for p in PRIMES]
rng = Random(8102026)
positive_gaps = 0
prime_power_gaps = 0

for case in range(80):
    exponents = [rng.randint(1, 3) for _ in PRIMES]
    if case % 2 == 0:
        rows = [tuple(bit * e for bit, e in zip(bits, exponents))
                for bits in product(range(2), repeat=3)]
    else:
        candidates = list(product(*(range(e + 1) for e in exponents)))
        rows = rng.sample(candidates, 8)
    points = [gprod(mul(gp(p, a), gp(conj(p), e - a))
                    for p, e, a in zip(PRIMES, exponents, row))
              for row in rows]
    order = sorted(range(8), key=lambda i: atan2(points[i][1], points[i][0]))
    points = [points[i] for i in order]
    rows = [rows[i] for i in order]
    assert len(set(points)) == 8
    assert len({norm(z) for z in points}) == 1

    chords, cores, residues = {}, {}, {}
    for i, j in combinations(range(8), 2):
        chord = divide(add(points[i], tuple(-t for t in points[j])), (0, 2))
        core = gprod(mul(gp(p, min(rows[i][k], rows[j][k])),
                         gp(conj(p), e - max(rows[i][k], rows[j][k])))
                     for k, (p, e) in enumerate(zip(PRIMES, exponents)))
        residue = divide(chord, core)
        assert residue[1] == 0 and residue[0] != 0
        chords[i, j], cores[i, j], residues[i, j] = chord, core, abs(residue[0])

    matching_edges = [[(i, 4 + s[i]) for i in range(4)] for s in PERMS]
    matching_products = [gprod(chords[edge] for edge in edges)
                         for edges in matching_edges]
    lcd = glcm(matching_products)
    ratio = divide(lcd, matching_products[0])
    if ratio[1] == 0:
        lcd = mul(lcd, (1 if ratio[0] > 0 else -1, 0))
    else:
        assert ratio[0] == 0
        lcd = mul(lcd, (0, -1 if ratio[1] > 0 else 1))
    summands = [divide(lcd, value) for value in matching_products]
    assert all(a[1] == 0 and a[0] > 0 for a in summands)
    summands = [a[0] for a in summands]
    content = 0
    for a in summands:
        content = gcd(content, a)
    assert content == 1
    permanent = sum(summands)
    signed_det = sum(s * a for s, a in zip(SIGNS, summands))
    determinant = abs(signed_det)
    assert determinant > 0
    # Borchardt after multiplication by the square of the actual Gaussian LCD.
    assert sum(s * a * a for s, a in zip(SIGNS, summands)) == signed_det * permanent

    lcd_core = (1, 0)
    forced_content = 1
    separated_factor = 1
    for k, (p, pn, e) in enumerate(zip(PRIMES, PNORMS, exponents)):
        first = max(sum(min(rows[i][k], rows[j][k]) for i, j in edges)
                    for edges in matching_edges)
        second = max(sum(e - max(rows[i][k], rows[j][k]) for i, j in edges)
                     for edges in matching_edges)
        layer_first = layer_second = content_order = 0
        for level in range(1, e + 1):
            a = sum(rows[i][k] >= level for i in range(4))
            b = sum(rows[i][k] >= level for i in range(4, 8))
            layer_first += min(a, b)
            layer_second += min(4 - a, 4 - b)
            content_order += comb(abs(a - b), 2)
        assert (first, second) == (layer_first, layer_second)
        lcd_core = mul(lcd_core, mul(gp(p, first), gp(conj(p), second)))
        forced_content *= pn ** content_order
        gap = max(0, min(rows[j][k] for j in range(4, 8))
                  - max(rows[i][k] for i in range(4)),
                  min(rows[i][k] for i in range(4))
                  - max(rows[j][k] for j in range(4, 8)))
        separated_factor *= pn ** gap
        positive_gaps += gap > 0
        prime_power_gaps += gap > 1
    quotient = divide(lcd, lcd_core)
    assert quotient[0] == 0 or quotient[1] == 0
    correction = abs(quotient[0] or quotient[1])
    cross = prod(residues[i, j] for i in range(4) for j in range(4, 8))
    within = prod(residues[i, j] for group in (range(4), range(4, 8))
                  for i, j in combinations(group, 2))
    assert cross % correction == 0
    assert determinant * cross == forced_content * correction * within
    assert all((a - summands[0]) % separated_factor == 0 for a in summands)
    assert gcd(permanent, separated_factor) == 1
    assert permanent >= separated_factor + 24
    assert (determinant // gcd(determinant, permanent)) % separated_factor ** 6 == 0

assert positive_gaps and prime_power_gaps
assert sum(comb(4, a) * comb(4, b) * abs(a - b)
           for a, b in product(range(5), repeat=2)) == 280
assert sum(comb(4, a) * comb(4, b) * comb(abs(a - b), 2)
           for a, b in product(range(5), repeat=2)) == 116
print("80 actual equal-norm Gaussian eight-tuples checked exactly.")
print("Cauchy residue quotient, Borchardt identity, positive primitive summands pass.")
print(f"Separating-prime cases: {positive_gaps}; gaps exceeding one: {prime_power_gaps}.")
print("Permanent coprimality, positive lower bound, sixth-power denominator all pass.")
print("These test tuples are not asserted to be endpoint clusters.")
