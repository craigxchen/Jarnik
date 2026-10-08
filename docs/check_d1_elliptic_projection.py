"""Finite checks for boundary pullbacks and the F_2^3 involution cap."""

from itertools import combinations, product


def canonical_boundary(n, subset):
    subset = frozenset(subset)
    complement = frozenset(range(n)) - subset
    candidates = (subset, complement)
    return min(tuple(sorted(s)) for s in candidates)


def canonical_boundary_on(universe, subset):
    universe = frozenset(universe)
    subset = frozenset(subset)
    complement = universe - subset
    return min(tuple(sorted(s)) for s in (subset, complement))


def boundaries(n):
    return boundary_set(range(n))


def boundary_set(universe):
    universe = tuple(universe)
    return {
        canonical_boundary_on(universe, S)
        for size in range(2, len(universe) - 1)
        for S in combinations(universe, size)
    }


def forget_one_pullback(n, omitted):
    retained = [i for i in range(n) if i != omitted]
    target = boundary_set(retained)
    pullbacks = {}
    for A in target:
        A = set(A)
        source_terms = {
            canonical_boundary(n, A),
            canonical_boundary(n, A | {omitted}),
        }
        pullbacks[tuple(sorted(A))] = source_terms
    return retained, target, pullbacks


def forget_two_pullback(n, omitted):
    retained = [i for i in range(n) if i not in omitted]
    target = boundary_set(retained)
    pullbacks = {}
    for A in target:
        A = set(A)
        source_terms = {
            canonical_boundary(n, A | set(U))
            for bits in product((0, 1), repeat=2)
            for U in ({omitted[j] for j, bit in enumerate(bits) if bit},)
        }
        pullbacks[tuple(sorted(A))] = source_terms
    return retained, target, pullbacks


def check_boundary_pullbacks():
    assert len(boundaries(6)) == 25
    assert len(boundaries(5)) == 10
    assert len(boundaries(4)) == 3

    # Forgetting one mark: each M_0,5 boundary has two lifts. Exactly the
    # five source boundaries {omitted,k} contract after forgetting omitted.
    for omitted in range(6):
        _, target, pullbacks = forget_one_pullback(6, omitted)
        assert len(target) == 10
        assert all(len(terms) == 2 for terms in pullbacks.values())
        used = set().union(*pullbacks.values())
        assert len(used) == 20
        contracted = {
            canonical_boundary(6, {omitted, k})
            for k in range(6) if k != omitted
        }
        assert len(contracted) == 5
        assert used == boundaries(6) - contracted

    # Forgetting two marks: each of the three M_0,4 boundary points has four
    # lifts, one for every placement of the two omitted labels.
    for omitted in combinations(range(6), 2):
        _, target, pullbacks = forget_two_pullback(6, omitted)
        assert len(target) == 3
        assert all(len(terms) == 4 for terms in pullbacks.values())
        assert len(set().union(*pullbacks.values())) == 12


def f2_add(x, y):
    return x ^ y


def is_cap(points):
    points = set(points)
    return all(f2_add(f2_add(x, y), z) != 0
               for x, y, z in combinations(points, 3))


def check_f2_cap():
    nonzero = tuple(range(1, 8))
    caps = [set(S) for size in range(8)
            for S in combinations(nonzero, size) if is_cap(S)]
    assert max(map(len, caps)) == 4
    assert any(S == {1, 2, 4, 7} for S in caps)
    assert all(len(S) <= 4 for S in caps)


if __name__ == "__main__":
    check_boundary_pullbacks()
    check_f2_cap()
    print("M_0,6 boundary count 25; one-mark pullbacks have 2 terms and two-mark pullbacks 4 terms.")
    print("F_2^3 cap enumeration: maximum 4 distinct nonzero involutions with no dependent triple.")
