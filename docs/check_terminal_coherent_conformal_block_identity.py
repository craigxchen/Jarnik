"""Finite coefficient audit of the terminal highest-root identity.

This checks explicit determinant products at d=1,2,4 and several
integer z-tuples. It does not prove the all-d identity or a height bound.
"""

from collections import defaultdict
from itertools import permutations
from math import factorial


def permutation_sign(perm):
    inversions = sum(perm[i] > perm[j]
                     for i in range(len(perm)) for j in range(i + 1, len(perm)))
    return -1 if inversions % 2 else 1


def triangle_product_terms(triples):
    """Expand a product of disjoint row determinants by coordinate words."""
    m = 3 * len(triples)
    terms = {(-1,) * m: 1}
    for triple in triples:
        updated = defaultdict(int)
        for word, coefficient in terms.items():
            for perm in permutations((0, 1, 2)):
                new_word = list(word)
                for label, coordinate in zip(triple, perm):
                    new_word[label] = coordinate
                updated[tuple(new_word)] += coefficient * permutation_sign(perm)
        terms = updated
    assert all(-1 not in word for word in terms)
    return terms


def check(triples, z):
    d = len(triples)
    m = 3 * d
    assert len(z) == m
    terms = triangle_product_terms(triples)
    coherent = defaultdict(int)
    derivative = defaultdict(int)
    for word, coefficient in terms.items():
        c_set = tuple(i for i, coord in enumerate(word) if coord == 2)
        b_set = tuple(i for i, coord in enumerate(word) if coord == 1)
        assert len(c_set) == len(b_set) == d
        z_on_b = 1
        for i in b_set:
            z_on_b *= z[i]
        coherent[c_set] += coefficient * z_on_b
        z_on_c = 1
        for i in c_set:
            z_on_c *= z[i]
        derivative[b_set] += factorial(d) * coefficient * z_on_c
        swapped = tuple(2 if c == 1 else 1 if c == 2 else 0 for c in word)
        assert terms[swapped] == (-1)**d * coefficient
    for subset in set(coherent) | set(derivative):
        assert derivative[subset] == (-1)**d * factorial(d) * coherent[subset]


def main():
    fixtures = [
        ((0, 1, 2),),
        ((0, 1, 2), (3, 4, 5)),
        ((0, 3, 5), (1, 2, 4)),
        ((0, 4, 8), (1, 5, 9), (2, 6, 10), (3, 7, 11)),
    ]
    for triples in fixtures:
        m = 3 * len(triples)
        check(triples, tuple(range(m)))
        check(triples, tuple((i * i + 3 * i + 1) % 7 for i in range(m)))
    print("Highest-root/coherent coefficient identity verified for d=1,2,4 fixtures.")


if __name__ == "__main__":
    main()
