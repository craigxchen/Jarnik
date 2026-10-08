"""Exact finite checks for the invariant-relation lattice note.

Checks dimensions, balanced ranks and Gram counts, attained core orders,
unimodular kernel coordinates, and the CRT relation-lattice identities.
The CRT vector is not asserted to be an actual invariant evaluation.
"""

from fractions import Fraction
from itertools import combinations, product
from math import comb, factorial, gcd, prod


def rank(a):
    a = [[Fraction(x) for x in row] for row in a]
    r = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        scale = a[r][j]
        a[r] = [x / scale for x in a[r]]
        for i in range(r + 1, len(a)):
            if a[i][j]:
                scale = a[i][j]
                a[i] = [x - scale * y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def matchings(labels):
    if not labels:
        yield []
        return
    a = labels[0]
    for j, b in enumerate(labels[1:], 1):
        for rest in matchings(labels[1:j] + labels[j + 1:]):
            yield [(a, b)] + rest


def basis(m):
    ans = []
    for top_counts in product(range(3), repeat=m):
        if sum(top_counts) != m:
            continue
        top = [i for i, c in enumerate(top_counts) for _ in range(c)]
        bottom = [i for i, c in enumerate(top_counts) for _ in range(2 - c)]
        if all(a < b for a, b in zip(top, bottom)):
            ans.append(list(zip(top, bottom)))
    return ans


def dimension(ds):
    coefficients = [1]
    for d in ds:
        nxt = [0] * (len(coefficients) + d)
        for i, c in enumerate(coefficients):
            for j in range(d + 1):
                nxt[i + j] += c
        coefficients = nxt
    total = sum(ds)
    return coefficients[total // 2] - coefficients[total // 2 - 1]


def evaluate(graph, inside):
    value = 1
    for a, b in graph:
        if (a in inside) == (b in inside):
            return 0
        value *= 1 if a in inside else -1
    return value


def internal(graph, inside):
    return sum(a in inside and b in inside for a, b in graph)


def triangle(v):
    a, b, c = sorted(v)
    return [(a, b), (b, c), (a, c)]


def kernel_graphs(m):
    labels = set(range(m))
    for six in combinations(range(m), 6):
        six = set(six)
        first = min(six)
        for other in combinations(sorted(six - {first}), 2):
            t = {first, *other}
            for rest in matchings(sorted(labels - six)):
                yield triangle(t) + triangle(six - t) + rest + rest


def egcd(a, b):
    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if old_r < 0:
        old_r, old_s, old_t = -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def column_kernel_coordinates(c):
    a = [row[:] for row in c]
    rows, cols = len(a), len(a[0])
    transform = [[int(i == j) for j in range(cols)] for i in range(cols)]
    for i in range(rows):
        pivot = next(j for j in range(i, cols) if a[i][j])
        for matrix in (a, transform):
            for row in matrix:
                row[i], row[pivot] = row[pivot], row[i]
        for j in range(i + 1, cols):
            if not a[i][j]:
                continue
            x, y = a[i][i], a[i][j]
            g, u, v = egcd(x, y)
            for matrix in (a, transform):
                for row in matrix:
                    old_i, old_j = row[i], row[j]
                    row[i] = u * old_i + v * old_j
                    row[j] = -(y // g) * old_i + (x // g) * old_j
    return a, transform


def determinant(a):
    a = [[Fraction(x) for x in row] for row in a]
    result = Fraction(1)
    for j in range(len(a)):
        pivot = next((i for i in range(j, len(a)) if a[i][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            result = -result
        x = a[j][j]
        result *= x
        for i in range(j + 1, len(a)):
            scale = a[i][j] / x
            a[i] = [b - scale * c for b, c in zip(a[i], a[j])]
    assert result.denominator == 1
    return int(result)


def inverse(a):
    n = len(a)
    a = [[Fraction(x) for x in row] + [Fraction(i == j) for j in range(n)]
         for i, row in enumerate(a)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        scale = a[j][j]
        a[j] = [x / scale for x in a[j]]
        for i in range(n):
            if i != j and a[i][j]:
                scale = a[i][j]
                a[i] = [x - scale * y for x, y in zip(a[i], a[j])]
    return [row[n:] for row in a]


saved_c = None
for m in (4, 6, 8):
    q = m // 2
    graphs = basis(m)
    cuts = [set(c) for c in combinations(range(m), q) if 0 in c]
    c = [[evaluate(g, cut) for g in graphs] for cut in cuts]
    assert len(graphs) == dimension([2] * m)
    assert rank(c) == len(cuts) == comb(m, q) // 2
    assert all(evaluate(g, cut) == evaluate(g, set(range(m)) - cut)
               for g in graphs for cut in cuts)
    ms = list(matchings(list(range(m))))
    incidence = [[int(all((a in cut) != (b in cut) for a, b in mat))
                  for mat in ms] for cut in cuts]
    gram = [[sum(x * y for x, y in zip(a, b)) for b in incidence]
            for a in incidence]
    assert all(gram[i][j] == factorial(len(s & t)) * factorial(q - len(s & t))
               for i, s in enumerate(cuts) for j, t in enumerate(cuts))
    assert rank(gram) == len(cuts)
    if m >= 6:
        kg = list(kernel_graphs(m))
        assert all(evaluate(g, cut) == 0 for g in kg for cut in cuts)
        for mask in range(1 << m):
            inside = {i for i in range(m) if mask >> i & 1}
            g0 = max(0, 2 * len(inside) - m)
            assert min(2 * internal(mat, inside) for mat in ms) == g0
            assert min(internal(g, inside) for g in kg) == g0 + (len(inside) == q)
    if m == 6:
        saved_c = c

assert saved_c is not None
cnew, transform = column_kernel_coordinates(saved_c)
s, r = len(cnew), len(cnew[0])
assert abs(determinant(transform)) == 1
assert all(cnew[i][j] == sum(saved_c[i][k] * transform[k][j] for k in range(r))
           for i in range(s) for j in range(r))
assert all(x == 0 for row in cnew for x in row[s:])
d = [row[:s] for row in cnew]
assert all(gcd(*row) == 1 for row in d)
det_d = determinant(d)
assert det_d
adj = [[x * det_d for x in row] for row in inverse(d)]
assert all(x.denominator == 1 for row in adj for x in row)
adj = [[int(x) for x in row] for row in adj]
moduli = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
modulus = prod(moduli)
u = [pow(modulus // p, -1, p) for p in moduli]
idempotents = [(modulus // p) * v for p, v in zip(moduli, u)]
n = 10**6
f = [sum(idempotents[j] * d[j][i] for j in range(s)) for i in range(s)]
f += [modulus * n, modulus] + [0] * (r - s - 2)
assert gcd(*f) == 1
assert gcd(*f[s:]) == modulus
for j, p in enumerate(moduli):
    assert all((f[i] - idempotents[j] * cnew[j][i]) % p == 0 for i in range(r))
    relation = [p * adj[i][j] for i in range(s)] + [0, -u[j] * det_d]
    relation += [0] * (r - s - 2)
    assert sum(a * b for a, b in zip(f, relation)) == 0
    assert any(sum(a * b for a, b in zip(row, relation)) for row in cnew)
    bound = p * max(abs(det_d), max(abs(x) for row in adj for x in row))
    assert max(map(abs, relation)) <= bound

print('Passed: exact invariant dimensions and balanced ranks for m=4,6,8.')
print('Passed: perfect-matching Gram counts, complementary coefficients, and both core minima.')
print('Passed: unimodular balanced-kernel coordinates and every CRT outside-relation identity.')
for m in (4, 6, 8, 10, 12, 14):
    r = dimension([2] * m)
    s = comb(m, m // 2) // 2
    a = m * 2**(m - 2) - m * s
    denominator = s if r > s else r - 1
    print(f'm={m}: r={r}, s={s}, A={a}, critical exponent={Fraction(a, denominator)}')
print('The CRT evaluation vector is an integer linear model, not a circle construction.')
