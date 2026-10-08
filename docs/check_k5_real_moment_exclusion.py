"""Exact finite certificates for real K5 first/third-moment exclusion.

No symbolic algebra or LP package is needed. The mathematical reduction to
four sign-order cones uses the six-root theorem cited in the companion note.
This verifies the exhaustive orders, complete linear parametrizations and
strictly positive polynomial certificates with integer arithmetic.
"""

from fractions import Fraction
from itertools import combinations, permutations

EDGES = tuple(combinations(range(5), 2))
PATTERNS = ((1, -1, -1, 1, 1, -1), (-1, 1, 1, -1, -1, 1),
            (1, 1, 1, -1, -1, 1), (-1, -1, -1, 1, 1, -1))
CONTRIBUTIONS = tuple(tuple(int(i in e) - int(j in e)
                            for i, j in EDGES) for e in EDGES)
# Signed edges, in order of increasing absolute value.
WORDS = (
    ((0, 1), (1, -1), (2, -1), (4, 1), (8, -1), (6, -1),
     (9, 1), (5, -1), (7, -1), (3, -1)),
    ((0, 1), (1, -1), (4, 1), (2, -1), (8, -1), (6, -1),
     (9, 1), (5, -1), (7, -1), (3, -1)),
    ((0, 1), (1, 1), (2, 1), (4, -1), (8, 1), (6, 1),
     (9, -1), (5, 1), (7, 1), (3, 1)),
    ((0, 1), (1, 1), (4, -1), (2, 1), (8, 1), (6, 1),
     (9, -1), (5, 1), (7, 1), (3, 1)),
)
# Each entry is a six-vector list of columns: g=sum_j lambda_j R_j.
COMMON_RAYS = (
    (0, 0, 0, 0, 0, 1, 0, 0, 1, 0),
    (0, 0, 1, 0, 1, 0, 1, 0, 0, 0),
    (0, 1, 0, 0, 1, 1, 0, 1, 0, 0),
    (0, 0, 0, 0, 1, 0, 1, 1, 0, 1),
)
RAYS = (
    ((1, 0, 0, 0, 0, 2, 0, 0, 0, 0),
     (0, 0, 0, 1, 0, 0, 1, 0, 0, 0), COMMON_RAYS[0],
     *COMMON_RAYS[1:]),
    ((0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
     (1, 0, 0, 0, 0, 2, 0, 0, 0, 0), COMMON_RAYS[0],
     *COMMON_RAYS[1:]),
    ((0, 0, 0, 1, 0, 0, 1, 0, 0, 0), COMMON_RAYS[0],
     (1, 0, 0, 0, 2, 0, 0, 2, 0, 0), *COMMON_RAYS[1:]),
    ((0, 0, 0, 1, 0, 0, 0, 0, 0, 0), COMMON_RAYS[0],
     (1, 0, 0, 0, 2, 0, 0, 2, 0, 0), *COMMON_RAYS[1:]),
)
# Constant (degree zero) or linear multipliers of F_1,...,F_4.
CERTIFICATES = (
    ((0,), (0,), (-1,), (3,)),
    ((0, 0, 0, 36, 36, 36), (0, 0, 36, 0, 346, 36),
     (-18, -36, 27, -364, -1266, -216),
     (18, 36, -27, 18, 920, 216)),
    ((-1,), (-1,), (1,), (-1,)),
    ((0, -2340, -11172, -3756, -5586, -312),
     (0, 2028, 4308, 3444, 5598, -312),
     (156, -1494, 208, 3276, 2871, 1872),
     (-156, 3834, 6032, -156, 2079, -1872)),
)
ZERO = (0,) * 6
UNITS = tuple(tuple(int(i == j) for i in range(6)) for j in range(6))


def admissible_words():
    leaves = []
    visits = [0] * 11

    def visit(word, used, lengths, masks):
        visits[len(word)] += 1
        if len(word) == 10:
            leaves.append(tuple(word))
            return
        for e in range(10):
            if used & (1 << e):
                continue
            for sign in (1, -1):
                if not word and (e != 0 or sign != 1):
                    continue
                ls, ms = list(lengths), list(masks)
                for p, contribution in enumerate(CONTRIBUTIONS[e]):
                    if not contribution:
                        continue
                    ms[p] = sum(1 << k for k in range(4)
                                if masks[p] & (1 << k)
                                and PATTERNS[k][ls[p]] == sign * contribution)
                    if not ms[p]:
                        break
                    ls[p] += 1
                else:
                    visit(word + [(e, sign)], used | (1 << e), ls, ms)
    visit([], 0, [0] * 10, [15] * 10)
    assert visits == [1, 1, 18, 108, 192, 216, 216, 192, 240, 48, 48]
    return leaves


def canonical(word):
    return min(tuple((EDGES.index(tuple(sorted(p[v] for v in EDGES[e]))), sign)
                     for e, sign in word)
               for p in permutations(range(5)) if set(p[:2]) == {0, 1})


def rank(rows):
    a = [list(map(Fraction, row)) for row in rows]
    r = 0
    for col in range(len(a[0])):
        pivot = next((j for j in range(r, len(a)) if a[j][col]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        scale = a[r][col]
        a[r] = [x / scale for x in a[r]]
        for j in range(len(a)):
            if j != r:
                scale = a[j][col]
                a[j] = [x - scale * y for x, y in zip(a[j], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def add(a, b, factor=1):
    out = dict(a)
    for m, x in b.items():
        out[m] = out.get(m, 0) + factor * x
        if out[m] == 0:
            del out[m]
    return out


def multiply(a, b):
    out = {}
    for u, x in a.items():
        for v, y in b.items():
            m = tuple(i + j for i, j in zip(u, v))
            out[m] = out.get(m, 0) + x * y
    return {m: x for m, x in out.items() if x}


def linear(values):
    return {m: x for m, x in zip(UNITS, values) if x}


def main():
    words = admissible_words()
    assert len(words) == len(set(words)) == 48
    assert sorted({canonical(word) for word in words}) == list(WORDS)
    assert all(sum(canonical(word) == rep for word in words) == 12 for rep in WORDS)
    print('PASS: exhaustive successor-pattern search, 48 orders in four orbits')
    for number, (word, rays, certificate) in enumerate(zip(WORDS, RAYS, CERTIFICATES)):
        v = [[sign * (int(vertex in EDGES[e]) - int(0 in EDGES[e]))
              for e, sign in word] for vertex in range(1, 5)]
        b = [[sum(row[j:]) for j in range(10)] for row in v]
        assert rank(b) == 4 and rank(rays) == 6
        assert all(sum(x * y for x, y in zip(row, ray)) == 0
                   for row in b for ray in rays)
        # Six literal gap coordinates recover lambda_j, so positive gaps
        # imply every lambda_j>0. This proves completeness without an LP.
        inverse = []
        for j in range(6):
            i = next(i for i in range(10)
                     if tuple(ray[i] for ray in rays) == UNITS[j])
            inverse.append(i)
        g = [linear([ray[i] for ray in rays]) for i in range(10)]
        r, partial = [], {}
        for gap in g:
            partial = add(partial, gap)
            r.append(partial)
        cubic = [multiply(multiply(p, p), p) for p in r]
        fs = []
        for row in v:
            f = {}
            for weight, p in zip(row, cubic):
                f = add(f, p, weight)
            fs.append(f)
        positive = {}
        for weights, f in zip(certificate, fs):
            multiplier = {ZERO: weights[0]} if len(weights) == 1 else linear(weights)
            positive = add(positive, multiply(multiplier, f))
        assert positive and all(x > 0 for x in positive.values())
        degree = 3 if len(certificate[0]) == 1 else 4
        assert all(sum(m) == degree for m in positive)
        print(f'PASS: orbit {number}, lambda gaps {inverse}, degree {degree}, '
              f'{len(positive)} positive coefficients, minimum {min(positive.values())}')
    print('PASS: no real nonzero distinct-magnitude K5 first/third-moment solution')


if __name__ == '__main__':
    main()
