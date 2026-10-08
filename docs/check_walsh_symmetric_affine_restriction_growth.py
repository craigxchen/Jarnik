#!/usr/bin/env python3
"""Exact binary restriction and physical-prime dictionary checks."""

import random
from collections import Counter


def dot(x, y):
    return bin(x & y).count("1") & 1


def linear(cols, x):
    out = 0
    for j, value in enumerate(cols):
        if (x >> j) & 1:
            out ^= value
    return out


def rank(vectors):
    pivots = {}
    for value in vectors:
        x = value
        while x:
            j = x.bit_length() - 1
            if j not in pivots:
                pivots[j] = x
                break
            x ^= pivots[j]
    return len(pivots)


def form(cols, x, y):
    return dot(x, linear(cols, y))


def restrict_label(basis, label):
    return sum(dot(u, label) << j for j, u in enumerate(basis))


def gram_columns(cols, basis):
    return [restrict_label(basis, linear(cols, u)) for u in basis]


def nondegenerate_complement(cols, basis):
    """Split orthogonal lines and alternating planes, leaving the radical."""
    remaining = list(basis)
    chosen = []
    while remaining:
        index = next((i for i, u in enumerate(remaining)
                      if form(cols, u, u)), None)
        if index is not None:
            u = remaining.pop(index)
            chosen.append(u)
            remaining = [w ^ (u if form(cols, u, w) else 0)
                         for w in remaining]
            continue
        pair = next(((i, j) for i in range(len(remaining))
                     for j in range(i + 1, len(remaining))
                     if form(cols, remaining[i], remaining[j])), None)
        if pair is None:
            assert all(form(cols, u, w) == 0
                       for u in remaining for w in remaining)
            break
        i, j = pair
        u, v = remaining[i], remaining[j]
        chosen.extend((u, v))
        others = [w for z, w in enumerate(remaining) if z not in pair]
        remaining = [w ^ (u if form(cols, w, v) else 0)
                     ^ (v if form(cols, w, u) else 0) for w in others]
    assert rank(chosen + remaining) == len(basis)
    assert all(form(cols, u, w) == 0 for u in chosen for w in remaining)
    assert rank(gram_columns(cols, chosen)) == len(chosen)
    return chosen


def extraction(cols, offset, exceptions):
    t = len(cols)
    q = len(exceptions).bit_length()
    assert q <= t
    mask = (1 << q) - 1
    occupied = {x & mask for x in exceptions}
    x0 = next(x for x in range(1 << q) if x not in occupied)
    h_basis = [1 << j for j in range(q, t)]
    restricted_rank = rank(gram_columns(cols, h_basis))
    assert restricted_rank >= rank(cols) - 2 * q
    basis = nondegenerate_complement(cols, h_basis)
    assert len(basis) == restricted_rank
    n = len(basis)
    s_cols = gram_columns(cols, basis)
    d = restrict_label(basis, linear(cols, x0) ^ offset)
    shift = next(y for y in range(1 << n) if linear(s_cols, y) == d)
    selected = [x0 ^ linear(basis, y ^ shift) for y in range(1 << n)]
    assert len(set(selected)) == 1 << n
    assert not set(selected).intersection(exceptions)
    for y, x in enumerate(selected):
        assert restrict_label(basis, linear(cols, x) ^ offset) == linear(s_cols, y)
    assert (1 << t) <= (1 << (t - rank(cols) + 2 * q)) * len(selected)
    return basis, selected, s_cols


def random_symmetric(rng, t):
    cols = [0] * t
    for i in range(t):
        for j in range(i, t):
            if rng.randrange(2):
                cols[i] |= 1 << j
                cols[j] |= 1 << i
    return cols


def algebra_checks():
    rng = random.Random(20260915)
    count = 0
    for t in (4, 5, 6, 8, 10, 12):
        for _ in range(12):
            cols = random_symmetric(rng, t)
            offset = rng.randrange(1 << t)
            for size in (0, 1, 2, 3, 5):
                exceptions = set(rng.sample(range(1 << t), size))
                extraction(cols, offset, exceptions)
                count += 1
    # The explicit line plus alternating-plane conversion used in the proof.
    cols = [1, 4, 2]
    basis = [3, 5, 7]
    assert rank(basis) == 3
    assert gram_columns(cols, basis) == [1, 2, 4]
    return count


def physical_check(cols, offset, extra_exception):
    t, b = len(cols), 5
    size = 1 << t
    assignment = [linear(cols, x) ^ offset for x in range(size)]
    exceptions = {x for x, a in enumerate(assignment) if a == 0}
    for x in exceptions:
        assignment[x] = 1
    if extra_exception is not None:
        x = extra_exception
        exceptions.add(x)
        assignment[x] = 1 + (assignment[x] % (size - 1))
    assert max(Counter(assignment).values()) <= b
    basis, selected, s_cols = extraction(cols, offset, exceptions)
    n = len(basis)
    assert n >= 2
    designated = {a: [] for a in range(1, size)}
    for x, a in enumerate(assignment):
        designated[a].append(x)
    columns = []
    for a in range(1, size):
        rows = designated[a]
        for j in range(b):
            columns.append((a, rows[j] if j < len(rows) else None))
    multiplicity = Counter(restrict_label(basis, a) for a, _ in columns)
    for a in range(1, 1 << n):
        assert multiplicity[a] == b * (1 << (t - n))
    assert multiplicity[0] == b * ((1 << (t - n)) - 1)
    selected_lookup = {x: y for y, x in enumerate(selected)}
    for a, flipped_row in columns:
        restricted = restrict_label(basis, a)
        local_flip = selected_lookup.get(flipped_row)
        orientation = (-1) ** dot(a, selected[0])
        if local_flip is not None:
            assert restricted == linear(s_cols, local_flip)
        for y, x in enumerate(selected):
            actual_sign = (-1) ** (dot(a, x) + int(flipped_row == x))
            predicted = orientation * (-1) ** (dot(restricted, y)
                                                + int(local_flip == y))
            assert actual_sign == predicted
        if restricted == 0:
            # The zero coordinate is the only possible exceptional entry.
            assert local_flip in (None, 0)
            assert all((-1) ** (dot(a, x) + int(flipped_row == x))
                       == orientation for x in selected[1:])
    # A balanced coefficient vector avoiding zero cancels every constant column.
    c = [0] * len(selected)
    c[1], c[2] = 1, -1
    for a, flipped_row in columns:
        if restrict_label(basis, a) == 0:
            assert sum(c[y] * (-1) ** (dot(a, x) + int(flipped_row == x))
                       for y, x in enumerate(selected)) == 0
    return len(columns) * len(selected)


def main():
    algebra = algebra_checks()
    comparisons = 0
    fixtures = 0
    for t in (4, 5, 6):
        identity = [1 << j for j in range(t)]
        forms = [identity]
        if t % 2 == 0:
            forms.append([1 << (j ^ 1) for j in range(t)])
        for cols in forms:
            for offset in (0, 3):
                comparisons += physical_check(cols, offset, None)
                fixtures += 1
                if t >= 5:
                    comparisons += physical_check(cols, offset, 1)
                    fixtures += 1
    print("PASS: %d binary affine restrictions; %d literal physical models; "
          "%d exact entry comparisons and constant-column cancellations"
          % (algebra, fixtures, comparisons))


if __name__ == "__main__":
    main()
