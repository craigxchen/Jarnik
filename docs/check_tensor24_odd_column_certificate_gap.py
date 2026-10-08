"""Exact bounded audit of one odd-column positive H12 tensor H2 code.

This is a binary-code and integer-lattice fixture.  It does not assert an
odd-moment solution or a Gaussian endpoint realization.
"""

import json
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import gcd
from pathlib import Path


def paley12():
    q = 11
    squares = {i * i % q for i in range(1, q)}
    out = [[1] * 12]
    for i in range(q):
        out.append([1] + [
            -1 if i == j else (1 if (i - j) % q in squares else -1)
            for j in range(q)
        ])
    return out


def tensor24():
    a = paley12()
    b = ((1, 1), (1, -1))
    return [[a[i][j] * b[u][v] for j in range(12) for v in range(2)]
            for i in range(12) for u in range(2)]


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def inverse(a):
    n = len(a)
    x = [[Fraction(a[i][j]) for j in range(n)] +
         [Fraction(i == j) for j in range(n)] for i in range(n)]
    for k in range(n):
        pivot = next(i for i in range(k, n) if x[i][k])
        x[k], x[pivot] = x[pivot], x[k]
        scale = x[k][k]
        x[k] = [z / scale for z in x[k]]
        for i in range(n):
            if i != k and x[i][k]:
                scale = x[i][k]
                x[i] = [z - scale * w for z, w in zip(x[i], x[k])]
    return [row[n:] for row in x]


def reconstruction(s):
    # Put mu=(x_0,...,x_6,-sum x_i).  Then c=B^t x.
    b = [[s[i][j] - s[7][j] for j in range(22)] for i in range(7)]
    pivots = []
    rank = 0
    for j in range(22):
        trial = pivots + [j]
        a = [[Fraction(b[i][k]) for i in range(7)] for k in trial]
        new_rank = 0
        for column in range(7):
            pivot = next((i for i in range(new_rank, len(a))
                          if a[i][column]), None)
            if pivot is None:
                continue
            a[new_rank], a[pivot] = a[pivot], a[new_rank]
            scale = a[new_rank][column]
            a[new_rank] = [z / scale for z in a[new_rank]]
            for i in range(len(a)):
                if i != new_rank and a[i][column]:
                    scale = a[i][column]
                    a[i] = [z - scale * w
                            for z, w in zip(a[i], a[new_rank])]
            new_rank += 1
        if new_rank > rank:
            pivots.append(j)
            rank = new_rank
        if rank == 7:
            break
    assert rank == 7 and len(pivots) == 7
    m = [[b[i][j] for i in range(7)] for j in pivots]
    m_inv = inverse(m)
    c_map = [[sum(Fraction(b[i][j]) * m_inv[i][k] for i in range(7))
              for k in range(7)] for j in range(22)]
    denominator = 1
    for row in c_map:
        for value in row:
            denominator = (denominator * value.denominator //
                           gcd(denominator, value.denominator))
    numerators = [[int(value * denominator) for value in row]
                  for row in c_map]
    return pivots, denominator, numerators


def signed_l1_vectors(dimension, radius):
    vector = [0] * dimension

    def visit(i, left):
        if i == dimension:
            yield tuple(vector)
            return
        for value in range(-left, left + 1):
            vector[i] = value
            yield from visit(i + 1, left - abs(value))

    yield from visit(0, radius)


def main():
    fixture = json.loads(Path(__file__).with_name(
        "tensor24_odd_column_certificate_gap_fixture.json").read_text())
    h = tensor24()
    assert all(dot(x, y) == (24 if i == j else 0)
               for i, x in enumerate(h) for j, y in enumerate(h))

    # Complete normalized-row scan; no row signings are included.
    profiles = Counter()
    positive = odd_positive = 0
    for rows in combinations(range(24), 8):
        counts = [sum(h[i][j] == 1 for i in rows) for j in range(24)]
        balanced = counts.count(4)
        if balanced < 14:
            continue
        assert counts[1] == 4
        positive += 1
        profile = tuple(Counter(min(x, 8 - x) for x in counts)[i]
                        for i in range(5))
        profiles[profile] += 1
        if any(x % 2 for x in counts):
            odd_positive += 1
            assert profile == tuple(fixture["tensor_scan"]["odd_full_profile"])
    assert positive == fixture["tensor_scan"]["positive_row_sets"]
    assert odd_positive == fixture["tensor_scan"]["odd_column_positive_row_sets"]
    assert sum(profiles.values()) == positive

    rows = fixture["rows"]
    deleted = set(fixture["deleted_columns"])
    assert deleted == {0, 1}
    keep = [j for j in range(24) if j not in deleted]
    t = [h[i] for i in rows]
    s = [[row[j] for j in keep] for row in t]
    assert s == fixture["retained_matrix"]
    counts = [sum(row[j] == 1 for row in s) for j in range(22)]
    profile = tuple(Counter(min(x, 8 - x) for x in counts)[i]
                    for i in range(1, 5))
    assert profile == tuple(fixture["retained_minority_profile"])
    n1, n2, n3, n4 = profile
    defect = 9 * n1 + 4 * n2 + n3
    assert defect == 32 and defect + n4 - 44 == 3
    assert min(sum(a != b for a, b in zip(x, y))
               for x, y in combinations(s, 2)) == 11

    # The deleted balanced column is exactly the parity column.  Restoring it
    # and the constant column gives the orthogonal padded matrix, so K=0.
    parity = [(-1) ** sum(x == -1 for x in row) for row in s]
    assert parity == [row[1] for row in t]
    padded = [row + [p, 1] for row, p in zip(s, parity)]
    assert all(dot(x, y) == (24 if i == j else 0)
               for i, x in enumerate(padded) for j, y in enumerate(padded))

    pivots, denominator, c_map = reconstruction(s)
    assert pivots == fixture["pivot_columns"]
    assert denominator == fixture["pivot_reconstruction_denominator"]
    checked = short = 0
    for y in signed_l1_vectors(7, 10):
        checked += 1
        c = []
        for row in c_map:
            numerator = dot(row, y)
            if numerator % denominator:
                break
            c.append(numerator // denominator)
        else:
            if any(c) and sum(map(abs, c)) <= 10:
                short += 1
    assert checked == fixture["pivot_vectors_with_l1_at_most_10"]
    assert short == fixture["short_certificate_count"] == 0

    witness = fixture["sharp_witness"]
    lam = witness["lambda"]
    numerators = [dot(lam, column) for column in zip(*s)]
    assert witness["q"] > 0 and all(n % witness["q"] == 0 for n in numerators)
    c = [n // witness["q"] for n in numerators]
    assert sum(lam) == 0 and c == witness["c"]
    assert sum(map(abs, c)) == witness["l1"] == 11
    print("PASS: 735,471 tensor row sets; 1,320 odd-column positive profiles;")
    print("433,905 exact pivot vectors; certificate lattice minimum is 11.")


if __name__ == "__main__":
    main()
