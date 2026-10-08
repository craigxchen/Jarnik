#!/usr/bin/env python3
"""Exact F_2 and physical-column checks for the affine-plane certificate."""

import random


def dot(x, y):
    return bin(x & y).count("1") & 1


def linear(cols, x):
    out = 0
    while x:
        low = x & -x
        out ^= cols[low.bit_length() - 1]
        x ^= low
    return out


def rank(vectors):
    basis = {}
    for value in vectors:
        x = value
        while x:
            bit = x.bit_length() - 1
            if bit not in basis:
                basis[bit] = x
                break
            x ^= basis[bit]
    return len(basis)


def quadratic(cols, x):
    return dot(x, linear(cols, x))


def isotropic_direction(cols, t):
    size = 1 << t
    basis = []
    span = {0}
    for _ in range(t // 3):
        for v in range(1, size):
            if v in span or quadratic(cols, v):
                continue
            lv = linear(cols, v)
            if all(dot(u, lv) == dot(v, linear(cols, u)) == 0
                   for u in basis):
                basis.append(v)
                span |= {x ^ v for x in list(span)}
                break
        else:
            raise AssertionError("the isotropic extension was not found")
    direction = tuple(sorted(span))
    assert len(direction) == 1 << (t // 3)
    assert all(dot(r, linear(cols, s)) == 0
               for r in direction for s in direction)
    return direction, basis


def aligned_cosets(cols, b, t, direction, basis):
    size = 1 << t
    reps = sorted({min(x ^ r for r in direction) for x in range(size)})
    assert len(reps) == size // len(direction)
    restrictions = {tuple(dot(linear(cols, x), u) for u in basis)
                    for x in range(size)}
    assert len(restrictions) >= 2
    good = []
    for x in reps:
        labels = [linear(cols, x ^ r) ^ b for r in direction]
        assert all(all(dot(label, s) == dot(labels[0], s)
                       for s in direction) for label in labels)
        if any(dot(labels[0], s) for s in direction):
            assert all(label != 0 for label in labels)
            good.append(x)
    cosets = size // len(direction)
    assert len(good) in (cosets, cosets - cosets // len(restrictions))
    assert len(good) >= size // (2 * len(direction))
    return good


def still_aligned(labels, t, direction, reps):
    good = []
    for x in reps:
        values = [labels[x ^ r] for r in direction]
        if any(dot(values[0], s) for s in direction) and all(
            all(dot(value, s) == dot(values[0], s)
                for s in direction) for value in values
        ):
            good.append(x)
    return good


def physical_character(labels, t, direction, rep, copies=5):
    size = 1 << t
    counts = {}
    flipped = {}
    for x, label in enumerate(labels):
        assert 0 < label < size
        slot = counts.get(label, 0)
        assert slot < copies
        flipped[copies * (label - 1) + slot] = x
        counts[label] = slot + 1
    restriction = labels[rep]
    assert any(dot(restriction, r) for r in direction)
    c = {rep ^ r: 1 - 2 * dot(restriction, r) for r in direction}
    assert sum(c.values()) == 0
    assert sum(abs(x) for x in c.values()) == len(direction)
    qualifying = 0
    total = 0
    for label in range(1, size):
        sum_unflipped = sum(cx * (1 - 2 * dot(label, x))
                            for x, cx in c.items())
        assert sum_unflipped % 2 == 0
        old = sum_unflipped // 2
        if abs(old) == len(direction) // 2:
            qualifying += 1
        else:
            assert old == 0
        for slot in range(copies):
            j = copies * (label - 1) + slot
            val = old
            if j in flipped:
                x = flipped[j]
                val -= c.get(x, 0) * (1 - 2 * dot(label, x))
            assert abs(val) <= len(direction) // 2
            total += abs(val)
    assert qualifying == size // len(direction)
    assert total == copies * size // 2 - len(direction)
    return total


def matrix_fixtures():
    rng = random.Random(7331)
    checks = 0
    physical = 0
    for t in (6, 7, 8, 9, 10, 12):
        size = 1 << t
        for deficient in (False, True):
            for trial in range(18 if t <= 8 else (5 if t <= 10 else 1)):
                cols = [1 << j for j in range(t)]
                if deficient:
                    cols[-1] = 0
                for _ in range(6 * t):
                    i, j = rng.sample(range(t), 2)
                    cols[i] ^= cols[j]
                assert rank(cols) == t - int(deficient)
                if deficient:
                    image = {linear(cols, x) for x in range(size)}
                    b = next(x for x in range(size) if x not in image)
                    labels = [linear(cols, x) ^ b for x in range(size)]
                    assert all(label != 0 for label in labels)
                    assert max(labels.count(label) for label in set(labels)) == 2
                else:
                    b = rng.randrange(size)
                    labels = [linear(cols, x) ^ b for x in range(size)]
                direction, basis = isotropic_direction(cols, t)
                good = aligned_cosets(cols, b, t, direction, basis)
                threshold = size // (2 * len(direction))
                assert len(good) >= threshold

                # Spoil exactly threshold-minus-one disjoint good planes.
                modified = labels[:]
                coordinate = next(1 << j for j in range(t)
                                  if dot(1 << j, basis[0]))
                for rep in good[: threshold - 1]:
                    modified[rep] ^= coordinate
                surviving = still_aligned(modified, t, direction, good)
                assert surviving
                checks += 1

                if trial == 0:
                    if not deficient:
                        zero = labels.index(0)
                        labels[zero] = 1
                    assert all(label != 0 for label in labels)
                    assert max(labels.count(label) for label in set(labels)) <= 5
                    patched_good = still_aligned(labels, t, direction, good)
                    assert patched_good
                    for copies in (5, 9):
                        physical_character(labels, t, direction,
                                           patched_good[0], copies)
                        physical += 1
    return checks, physical


def three_dimensional_quadratics():
    checked = 0
    for cols0 in range(8):
        for cols1 in range(8):
            for cols2 in range(8):
                cols = [cols0, cols1, cols2]
                assert any(quadratic(cols, u) == 0 for u in range(1, 8))
                checked += 1
    return checked


if __name__ == "__main__":
    q = three_dimensional_quadratics()
    matrices, physical = matrix_fixtures()
    print(f"PASS: {q} binary three-space forms, {matrices} affine matrices "
          f"and exception checks, {physical} literal physical-copy characters")
