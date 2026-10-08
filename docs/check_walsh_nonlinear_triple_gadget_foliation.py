#!/usr/bin/env python3
"""Bounded exact checks for the nonlinear three-bit aligned foliation."""

from collections import Counter
from itertools import product


F = (0, 1, 3, 6, 7, 4, 5, 2)


def parity(x):
    return bin(x).count("1") & 1


def dot(x, y):
    return parity(x & y)


def label(row, blocks):
    out = 0
    for i in range(blocks):
        out |= F[(row >> (3 * i)) & 7] << (3 * i)
    return out


def restriction(lab, directions):
    return sum(dot(lab, v) << i for i, v in enumerate(directions))


def fixture(m, selected_choices):
    blocks = 2 * m
    M = 1 << (3 * blocks)
    h = 1 << m
    labels = [label(x, blocks) for x in range(M)]
    assert sorted(labels) == list(range(M))
    repaired = labels[:]
    repaired[0] = 1
    occ = Counter(repaired)
    assert occ[0] == 0 and max(occ.values()) == 2
    assert sum(occ.values()) == M

    for choices in selected_choices:
        assert len(choices) == m
        directions = []
        for i, (w0, w1) in enumerate(choices):
            assert 1 <= w0 <= 7 and 1 <= w1 <= 7
            directions.append((w0 << (6 * i)) | (w1 << (6 * i + 3)))
        V = {
            sum(v for bit, v in enumerate(directions) if mask & (1 << bit))
            for mask in range(h)
        }
        assert len(V) == h
        original_counts = Counter(restriction(a, directions) for a in labels)
        assert all(original_counts[ell] == M // h for ell in range(h))

        seen = set()
        coset_counts = Counter()
        for x in range(M):
            if x in seen:
                continue
            coset = {x ^ v for v in V}
            assert len(coset) == h
            seen.update(coset)
            ell = restriction(labels[x], directions)
            assert all(restriction(labels[y], directions) == ell for y in coset)
            if ell:
                assert all(restriction(repaired[y], directions) == ell for y in coset)
            coset_counts[ell] += 1
        assert len(seen) == M
        assert all(coset_counts[ell] == M // (h * h) for ell in range(h))
        assert sum(coset_counts[ell] * h for ell in range(1, h)) == M - M // h

    return M, h, M // (h * h)


def main():
    assert sorted(F) == list(range(8))
    for u, w in product(range(8), repeat=2):
        assert dot(F[u] ^ F[u ^ w], w) == int(w != 0)
    for w in range(1, 8):
        derivative = Counter(F[u] ^ F[u ^ w] for u in range(8))
        assert sorted(derivative.values()) == [2, 2, 2, 2]
    max_affine_agreement = 0
    for shift in range(8):
        for columns in product(range(8), repeat=3):
            agreement = 0
            for u in range(8):
                affine = shift
                for j in range(3):
                    if u & (1 << j):
                        affine ^= columns[j]
                agreement += affine == F[u]
            max_affine_agreement = max(max_affine_agreement, agreement)
    assert max_affine_agreement == 4
    one = [(w0, w1) for w0, w1 in product(range(1, 8), repeat=2)]
    assert fixture(1, [[choice] for choice in one]) == (64, 2, 16)
    two = [
        [(1, 2), (3, 5)],
        [(7, 4), (2, 6)],
        [(5, 1), (6, 3)],
    ]
    assert fixture(2, two) == (4096, 4, 256)
    print("nonlinear gadget and m=1,2 aligned foliations: OK")


if __name__ == "__main__":
    main()
