"""Exact finite audits of ordered sampling; no asymptotic proof by testing."""

from fractions import Fraction as F
from itertools import combinations, product
from math import comb
from random import Random


def prod(xs):
    ans = 1
    for x in xs:
        ans *= x
    return ans


def fixture(order, m, seed):
    rng = Random(seed)
    rows = list(range(order))
    rng.shuffle(rows)
    flips = [rng.choice((-1, 1)) for _ in range(1, order)]
    # Delete the constant Walsh column; orient columns independently.
    words = [tuple(flips[j - 1] * (-1) ** bin(i & j).count("1")
                   for i in rows[:m]) for j in range(1, order)]
    assert all(sum(w[i] * w[j] for w in words) == -1
               for i, j in combinations(range(m), 2))
    return words


def audit(order, m, k, seed):
    words = fixture(order, m, seed)
    positions = [tuple(a) for r in range(k + 1)
                 for a in combinations(range(k), r)]
    samples = list(combinations(range(m), k))
    mu = {(): F(1)}
    for r in range(1, k + 1):
        for j in combinations(range(m), r):
            mu[j] = F(sum(prod(w[i] for i in j) for w in words), len(words))

    energies = {}
    for r in range(1, k + 1):
        row_energy = sum(mu[j] ** 2 for j in combinations(range(m), r))
        assert row_energy <= F(2 * comb(m, r), m - r + 1)
        ordered_energy = sum(mu[tuple(i[a] for a in aa)] ** 2
                             for i in samples for aa in positions if len(aa) == r)
        ordered_energy /= len(samples)
        assert ordered_energy == F(comb(m - r, k - r), comb(m, k)) * row_energy
        assert ordered_energy <= F(2 * comb(k, r), m - r + 1)
        energies[r] = ordered_energy

    cube = list(product((-1, 1), repeat=k))
    target = tuple(1 if j in (0, k - 1) else -1 for j in range(k))
    tests = [
        ({e: F(e == target) for e in cube}, False),
        ({e: F(e == target or e == tuple(-x for x in target)) for e in cube}, True),
        ({e: F(e[0] * e[1] - 2 * e[1] * e[-1]) for e in cube}, True),
        ({e: F(3 * e[0] + e[1] * e[-1] - 2 * prod(e)) for e in cube}, False),
    ]
    for values, symmetric in tests:
        hat = {a: sum(values[e] * prod(e[j] for j in a) for e in cube) / (2 ** k)
               for a in positions}
        assert all(sum(hat[a] * prod(e[j] for j in a) for a in positions) == values[e]
                   for e in cube)
        if symmetric:
            assert all(hat[a] == 0 for a in positions if len(a) % 2)
        sample_values = []
        absolute_degrees = {r: F(0) for r in range(1, k + 1)}
        for i in samples:
            degree = {r: sum(hat[a] * mu[tuple(i[j] for j in a)]
                             for a in positions if len(a) == r)
                      for r in range(1, k + 1)}
            direct = sum(values[tuple(w[j] for j in i)] for w in words) / len(words)
            assert direct == hat[()] + sum(degree.values())
            sample_values.append(direct)
            for r in degree:
                absolute_degrees[r] += abs(degree[r]) / len(samples)
        for r in absolute_degrees:
            norm2 = sum(hat[a] ** 2 for a in positions if len(a) == r)
            assert absolute_degrees[r] ** 2 <= norm2 * energies[r]

        if symmetric:
            # k almost equal consecutive blocks; each sample can meet all k.
            sizes = [m // k + (b < m % k) for b in range(k)]
            blocks = []
            for b, size in enumerate(sizes):
                blocks.extend([b] * size)
            collision = F(sum(len({blocks[j] for j in i}) < k for i in samples), len(samples))
            assert collision <= F(comb(k, 2) * (max(sizes) - 1), m - 1)
            for b, size in enumerate(sizes):
                bias2 = sum(F(sum(w[j] for j in range(m) if blocks[j] == b), size) ** 2
                            for w in words) / len(words)
                assert bias2 <= F(1, size)
            error = abs(sum(sample_values) / len(samples) - hat[()])
            l1 = sum(abs(hat[a]) for a in positions if a)
            osc = max(values.values()) - min(values.values())
            assert error <= l1 / min(sizes) + osc * collision
    return len(samples)


def main():
    total = sum(audit(*args) for args in [(8, 8, 5, 31), (16, 12, 5, 79), (16, 16, 4, 11)])
    # A fixed rank subset need not have uniform row-set multiplicities.
    counts = {i: 0 for i in range(6)}
    for sample in combinations(range(6), 3):
        counts[sample[0]] += 1
    assert len(set(counts.values())) > 1
    assert F(1, 16) / F(15, 16) == F(1, 15)
    print(f"PASS: {total} ordered samples, four tests each; exact energy and block bounds")
    print("PASS: fixed-rank bias warning and intrinsic five-point ratio 1/15")


if __name__ == "__main__":
    main()
