"""Exact density-flat extraction and residual bookkeeping checks."""

from fractions import Fraction
from itertools import combinations
from random import Random


def span(basis):
    result = {0}
    for b in basis:
        result |= {v ^ b for v in tuple(result)}
    return result


def complement(t, vectors):
    generated = set(vectors)
    basis = []
    for j in range(t):
        b = 1 << j
        if b not in generated:
            basis.append(b)
            generated |= {x ^ b for x in tuple(generated)}
    return span(basis)


def extract(ambient, selected, target):
    t = ambient.bit_length() - 1
    vectors = {0}
    rows = set(selected)
    steps = 0
    while len(vectors) < target:
        reps = sorted(complement(t, vectors) & rows)
        n, h = len(reps), len(vectors)
        assert n * h == len(rows)
        if n < 2:
            return None, steps
        directions = {}
        for x, y in combinations(reps, 2):
            d = x ^ y
            directions.setdefault(d, []).append((x, y))
        d = max(directions, key=lambda x: len(directions[x]))
        pairs = directions[d]
        assert 4 * len(pairs) * ambient >= n * n * h
        new_vectors = vectors | {v ^ d for v in vectors}
        new_rows = {x ^ v for x, y in pairs for v in new_vectors}
        assert new_rows <= rows
        assert len(new_rows) == len(pairs) * 2 * h
        assert 2 * ambient * len(new_rows) >= len(rows) ** 2
        vectors, rows = new_vectors, new_rows
        steps += 1
    x0 = next(iter(rows))
    flat = {x0 ^ v for v in vectors}
    assert len(flat) == target and flat <= selected
    return flat, steps


def density_fixtures():
    rng = Random(260916)
    count = merges = 0
    fixtures = [(4, set(a), 2) for k in range(1, 5)
                for a in combinations(range(4), k)]
    for t in (5, 7, 9, 10):
        m = 1 << t
        for j in range(8):
            rows = set(rng.sample(range(m), (j + 8) * m // 16))
            fixtures.append((m, rows, 4))
    for m, rows, target in fixtures:
        flat, steps = extract(m, rows, target)
        # The density condition is checked without floating roots.
        sufficient = len(rows) ** target >= (2 ** target) * target * m ** (target - 1)
        if sufficient:
            assert flat is not None
        count += 1
        merges += steps
    return count, merges


def bookkeeping_fixtures():
    rng = Random(917)
    count = 0
    for m in (8, 16, 32):
        for j in range(20):
            wa = [Fraction(rng.randrange(5, 40), rng.randrange(1, 5))
                  for a in range(1, m)]
            fa = [w * Fraction(rng.randrange(0, 5), 5) for w in wa]
            w0, f = sum(wa), sum(fa)
            d = Fraction(rng.randrange(0, 5), 3)
            am = -Fraction(rng.randrange(1, 10), 2)
            ell = (w0 + d - 2 * f + 2 * am) / m
            r = [w - 4 * z / m - ell for w, z in zip(wa, fa)]
            mu = [x - 2 * f / m for x in r]
            qstar = (w0 + d + 2 * am) / m
            assert all(w == qstar + x + 4 * z / m
                       for w, x, z in zip(wa, mu, fa))
            negative = sum(max(0, -x) for x in mu)
            assert sum(mu) == (w0 - 4 * f) / m - Fraction(m - 1, m) * (d + 2 * am)
            assert d <= (w0 - 4 * f) / (m - 1) + Fraction(m, m - 1) * negative - 2 * am
            mean = w0 / (m - 1)
            deficit = sum(max(0, mean - w) for w in wa)
            assert deficit <= negative + (m - 1) * max(0, mean - qstar)
            assert sum(abs(w - mean) for w in wa) == 2 * deficit
            assert (m - 1) * max(0, mean - qstar) <= w0 / m + max(0, -2 * am)
            count += 1
    return count


if __name__ == "__main__":
    fixtures, merges = density_fixtures()
    bookkeeping = bookkeeping_fixtures()
    print("PASS: %d density fixtures; %d exact affine-coset merges; "
          "%d rational residual/class-mass identities."
          % (fixtures, merges, bookkeeping))
