"""Exact finite checks for the arbitrary-vector-field aligned-flat extraction."""

from collections import Counter
from itertools import product
from math import floor, log2
from random import Random


def dot(x, y):
    return bin(x & y).count('1') & 1


def rank_binary(matrix):
    pivots = {}
    for row in matrix:
        bits = sum(value << j for j, value in enumerate(row))
        while bits:
            pivot = bits.bit_length() - 1
            if pivot in pivots:
                bits ^= pivots[pivot]
            else:
                pivots[pivot] = bits
                break
    return len(pivots)


def span(basis):
    vectors = [0]
    for b in basis:
        vectors += [v ^ b for v in vectors]
    return tuple(sorted(vectors))


def complement(t, vectors):
    basis = []
    generated = set(vectors)
    for j in range(t):
        b = 1 << j
        if b not in generated:
            basis.append(b)
            generated |= {v ^ b for v in tuple(generated)}
    assert len(generated) == 1 << t
    return span(basis)


def aligned(labels, vectors, x, restriction):
    return all(tuple(dot(labels[x ^ v], w) for w in vectors) == restriction
               for v in vectors)


def start(labels, t):
    m = 1 << t
    q = [[dot(labels[x], x ^ y) for y in range(m)] for x in range(m)]
    aa = [[1 ^ q[x][y] ^ q[y][x] for y in range(m)] for x in range(m)]
    assert all(sum(row) == m // 2 for row in q)
    rr = rank_binary(aa)
    assert rr <= 2 * t + 3
    count00 = sum(q[x][y] == q[y][x] == 0
                  for x in range(m) for y in range(m))
    count11 = sum(q[x][y] == q[y][x] == 1
                  for x in range(m) for y in range(m))
    assert count00 == count11
    assert sum(map(sum, aa)) * rr >= m * m
    assert 2 * rr * count11 >= m * m
    directions = {
        w: [x for x in range(m)
            if dot(labels[x], w) == dot(labels[x ^ w], w) == 1]
        for w in range(1, m)
    }
    assert sum(map(len, directions.values())) == count11
    w = max(directions, key=lambda v: len(directions[v]))
    rows = set(directions[w])
    assert len(rows) * 2 * (2 * t + 3) >= m
    vectors = (0, w)
    restriction = (0, 1)
    reps = [x for x in complement(t, vectors) if x in rows]
    assert len(reps) * 2 == len(rows)
    assert all(aligned(labels, vectors, x, restriction) for x in reps)
    return vectors, restriction, reps


def extend(labels, t, vectors, restriction, representatives):
    m, h, n = 1 << t, len(vectors), len(representatives)
    comp = complement(t, vectors)
    old_rows = {x ^ v for x in representatives for v in vectors}
    reps = [x for x in comp if x in old_rows]
    assert len(reps) == n
    assert all(aligned(labels, vectors, x, restriction) for x in reps)
    kk = [[int(all(dot(labels[x ^ v], x ^ y) == 0
                        and dot(labels[y ^ v], x ^ y) == 0 for v in vectors))
           for y in reps] for x in reps]
    assert all(kk[j][j] == 1 for j in range(n))
    assert all(kk[i][j] == kk[j][i] for i in range(n) for j in range(n))
    rr = rank_binary(kk)
    assert rr <= (t + 2) ** (2 * h)
    assert sum(map(sum, kk)) * rr >= n * n
    pairs = Counter()
    for i, x in enumerate(reps):
        for j, y in enumerate(reps):
            if i != j and kk[i][j]:
                pairs[x ^ y] += 1
    assert sum(pairs.values()) == sum(map(sum, kk)) - n
    assert all(w in comp and w not in vectors for w in pairs)
    assert len(pairs) <= m // h - 1
    if not pairs:
        assert n < 2 * rr
        return None
    w = max(pairs, key=pairs.get)
    new_vectors = tuple(sorted(vectors + tuple(v ^ w for v in vectors)))
    old_value = dict(zip(vectors, restriction))
    new_restriction = tuple(old_value[v if v in old_value else v ^ w]
                            for v in new_vectors)
    eligible = {x for i, x in enumerate(reps) for j, y in enumerate(reps)
                if x ^ y == w and kk[i][j]}
    new_rows = {x ^ v for x in eligible for v in vectors}
    new_reps = [x for x in complement(t, new_vectors) if x in new_rows]
    assert len(new_reps) * 2 == pairs[w]
    assert len(new_rows) == len(new_reps) * 2 * h
    assert any(new_restriction)
    assert all(aligned(labels, new_vectors, x, new_restriction)
               for x in new_reps)
    if n >= 2 * rr:
        # Stronger fixture version with actual rank in place of its bound.
        assert len(new_rows) * 2 * rr * m >= len(old_rows) ** 2
    return new_vectors, new_restriction, new_reps


def physical_height(labels, vectors, restriction, x0, b):
    m, h = len(labels), len(vectors)
    copy = {}
    used = Counter()
    for x, a in enumerate(labels):
        copy[x] = used[a]
        used[a] += 1
    assert max(used.values()) <= b
    cc = {x0 ^ v: (-1) ** restriction[j] for j, v in enumerate(vectors)}
    assert sum(cc.values()) == 0
    exponents = []
    for a in range(1, m):
        baseline = sum(c * (-1) ** dot(a, x) for x, c in cc.items()) // 2
        qualifies = tuple(dot(a, v) for v in vectors) == restriction
        assert abs(baseline) == (h // 2 if qualifies else 0)
        for j in range(b):
            correction = sum(c * (-1) ** dot(a, x) for x, c in cc.items()
                             if labels[x] == a and copy[x] == j)
            exponents.append(baseline - correction)
    assert sum(map(abs, exponents)) == b * m // 2 - h
    assert any(exponents)
    # Rational near-equal weights check only the exact error calculation;
    # they are not substituted for actual Gaussian prime logs in a theorem.
    r = len(exponents)
    scaled_weights = [r * r * 10 + j for j in range(r)]
    excess = (2 * sum(abs(v) * w for v, w in zip(exponents, scaled_weights))
              - sum(scaled_weights))
    assert excess < (b - 2 * h) * 10 * r * r + 2 * r * r


def main():
    count = 0
    for labels in product(range(1, 4), repeat=4):
        vv, ll, reps = start(labels, 2)
        extend(labels, 2, vv, ll, reps)
        count += 1
    assert count == 81
    rng = Random(260915)
    fixture = (28,21,11,30,24,23,13,3,12,7,25,4,31,1,22,3,
               19,18,10,9,2,16,17,27,26,6,15,5,14,8,20,29)
    vv, ll, reps = start(fixture, 5)
    assert extend(fixture, 5, vv, ll, reps) is None
    physical_height(fixture, vv, ll, reps[0], 5)
    successful_merges = 0
    for t in range(3, 8):
        m = 1 << t
        for case in range(12):
            labels = list(range(1, m)) + [rng.randrange(1, m)]
            rng.shuffle(labels)
            vv, ll, reps = start(labels, t)
            physical_height(labels, vv, ll, reps[0], 5)
            while True:
                result = extend(labels, t, vv, ll, reps)
                if result is None:
                    break
                successful_merges += 1
                vv, ll, reps = result
                physical_height(labels, vv, ll, reps[0], 5)
    assert successful_merges > 0
    # A strongly nonlinear field with a prepared aligned-plane family
    # exercises the quotient calculation after the starting stage.
    t = 7
    vv, ll = (0, 1, 2, 3), (0, 1, 0, 1)
    labels = [1 | (rng.randrange(1 << (t - 2)) << 2) for x in range(1 << t)]
    reps = list(complement(t, vv))
    assert all(aligned(labels, vv, x, ll) for x in reps)
    result = extend(labels, t, vv, ll, reps)
    if result is not None:
        nv, nl, nr = result
        physical_height(labels, nv, nl, nr[0], max(Counter(labels).values()) + 5)
    for t in list(range(4096, 20000, 37)) + [10**k for k in range(5, 16)]:
        lt = log2(t + 2)
        d = floor(log2(t / (8 * lt * lt)))
        h = 1 << d
        assert d >= 1
        assert t / (16 * lt * lt) <= h <= t / (8 * lt * lt)
        assert d + h * (d + 2) * lt + h / 2 < t
        assert h * d * lt + h / 2 <= t / 4
        previous = 2 * lt
        for i in range(1, d):
            previous = 2 * previous + 1 + (1 << (i + 1)) * lt
            closed = (1 << (i + 1)) * (i + 1) * lt + (1 << i) - 1
            assert abs(previous - closed) <= 1e-10 * max(1, closed)
    print(f"passed: 81 exhaustive fields; nonlinear fixtures; {successful_merges} merges; exact heights; recurrence cutoffs")


if __name__ == '__main__':
    main()
