"""Exact checks for the eight-pattern signed-witness parity gap.

Checks the Gaussian product table and finite divisor bounds, the rank-two
integer kernel, finite Walsh-code correlation counts, the independent Walsh
family, and the low-cost monomial exponent classification.
"""

from itertools import combinations, product
from math import gcd
from math import prod
import random


Gaussian = tuple[int, int]


def gmul(z: Gaussian, w: Gaussian) -> Gaussian:
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def gconj(z: Gaussian) -> Gaussian:
    return z[0], -z[1]


def gnorm(z: Gaussian) -> int:
    return z[0] * z[0] + z[1] * z[1]


def det(u: Gaussian, v: Gaussian) -> int:
    return u[0] * v[1] - v[0] * u[1]


def dot3(u: tuple[int, int, int], v: tuple[int, int, int]) -> int:
    return sum(x * y for x, y in zip(u, v))


def cross3(u: tuple[int, int, int], v: tuple[int, int, int]) -> tuple[int, int, int]:
    return (u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0])


SIGN_ROWS = (
    (1, 1, 1, 1, 1, 1, 1, 1),
    (1, 1, -1, -1, 1, 1, -1, -1),
    (1, -1, 1, -1, 1, -1, 1, -1),
    (1, -1, -1, 1, -1, 1, 1, -1),
)


def gaussian_product(factors) -> Gaussian:
    z = (1, 0)
    for factor in factors:
        z = gmul(z, factor)
    return z


def sign_products(groups: tuple[Gaussian, ...], corrections: tuple[Gaussian, ...]) -> tuple[Gaussian, ...]:
    assert len(groups) == 8 and len(corrections) == 4
    rows = []
    for i, signs in enumerate(SIGN_ROWS):
        factors = tuple(g if s == 1 else gconj(g) for g, s in zip(groups, signs))
        rows.append(gmul(corrections[i], gaussian_product(factors)))
    return tuple(rows)


def pair_g_norms(groups: tuple[Gaussian, ...]) -> dict[tuple[int, int], int]:
    return {(i, j): prod(gnorm(groups[g]) for g in range(8)
                         if SIGN_ROWS[i][g] == SIGN_ROWS[j][g])
            for i, j in combinations(range(4), 2)}


def expected_g_norms(groups: tuple[Gaussian, ...]) -> dict[tuple[int, int], int]:
    e = tuple(gnorm(z) for z in groups[:4])
    o = tuple(gnorm(z) for z in groups[4:])
    return {
        (0, 1): e[0] * e[1] * o[0] * o[1],
        (2, 3): e[0] * e[1] * o[2] * o[3],
        (0, 2): e[0] * e[2] * o[0] * o[2],
        (1, 3): e[0] * e[2] * o[1] * o[3],
        (0, 3): e[0] * e[3] * o[1] * o[2],
        (1, 2): e[0] * e[3] * o[0] * o[3],
    }


def fixture_check(groups: tuple[Gaussian, ...], corrections: tuple[Gaussian, ...],
                  require_cutoff: bool = True) -> tuple[bool, int, int]:
    norms = tuple(gnorm(z) for z in groups)
    assert all(n > 0 for n in norms)
    assert all(gcd(a, b) == 1 for a, b in combinations(norms, 2))
    # These are split-prime Gaussian generators or units; gcd coordinates one
    # certifies conjugate primitivity for the odd prime-norm cases.
    assert all(gcd(abs(z[0]), abs(z[1])) == 1 for z in groups)
    assert all(gnorm(k) > 0 for k in corrections)
    e = norms[:4]
    o = norms[4:]
    P = prod(norms)
    Ymax = 0
    M2 = max(gnorm(k) for k in corrections)
    points = sign_products(groups, corrections)
    y = tuple(z[1] for z in points)
    Ymax = max(map(abs, y))
    delta = {(i, j): det(points[i], points[j]) for i, j in combinations(range(4), 2)}
    Gij = pair_g_norms(groups)
    assert Gij == expected_g_norms(groups)
    assert all(delta[pair] % Gij[pair] == 0 for pair in Gij)

    # M<Q_* iff M^2 G_ij<P for every pair. The theorem forces all listed
    # coordinates and determinants nonzero under this strict cutoff.
    cutoff = all(M2 * g < P for g in Gij.values())
    if require_cutoff:
        assert cutoff
    if cutoff:
        assert all(v != 0 for v in y)
        assert all(v != 0 for v in delta.values())

    all_nondzero = all(v != 0 for v in y) and all(v != 0 for v in delta.values())
    if not all_nondzero:
        assert not cutoff
        return False, 0, Ymax

    t = {pair: delta[pair] // Gij[pair] for pair in Gij}
    assert all(v != 0 for v in t.values())
    T = max(map(abs, t.values()))
    Gmin = min(Gij.values())
    Omax = max(o)

    rows = (
        (y[2] * t[(0, 1)] * o[1], -y[1] * t[(0, 2)] * o[2], y[0] * t[(1, 2)] * o[3]),
        (y[3] * t[(0, 1)] * o[0], y[0] * t[(1, 3)] * o[3], -y[1] * t[(0, 3)] * o[2]),
        (y[0] * t[(2, 3)] * o[3], y[3] * t[(0, 2)] * o[0], -y[2] * t[(0, 3)] * o[1]),
        (y[1] * t[(2, 3)] * o[2], -y[2] * t[(1, 3)] * o[1], y[3] * t[(1, 2)] * o[0]),
    )
    kernel = (e[1], e[2], e[3])
    assert all(dot3(row, kernel) == 0 for row in rows)
    crosses = [cross3(rows[i], rows[j]) for i, j in combinations(range(4), 2)]
    nonzero = [v for v in crosses if v != (0, 0, 0)]
    assert nonzero
    for v in nonzero:
        assert v[0] * e[2] == v[1] * e[1]
        assert v[0] * e[3] == v[2] * e[1]
        assert v[0] % e[1] == 0
        scalar = v[0] // e[1]
        assert scalar != 0 and v == tuple(scalar * z for z in kernel)

    assert max(e[1:]) <= 2 * Ymax ** 2 * T ** 2 * Omax ** 2
    # Exact square form of T <= 2 M sqrt(P)Y/G.
    assert T * T * Gmin * Gmin <= 4 * M2 * P * Ymax * Ymax
    assert max(e[1:]) * Gmin * Gmin <= 8 * M2 * P * Ymax ** 4 * Omax ** 2
    return True, T, Ymax


def check_gaussian_fixtures() -> tuple[int, int, int, int]:
    prime_norm_groups = ((2, 1), (3, 2), (4, 1), (5, 2),
                         (6, 1), (5, 4), (7, 2), (6, 5))
    # Norms 5,13,17,29,37,41,53,61 are distinct split primes.
    assert tuple(map(gnorm, prime_norm_groups)) == (5, 13, 17, 29, 37, 41, 53, 61)
    correction_pool = ((1, 0), (1, 1), (1, -1), (2, 1), (2, -1), (0, 1), (2, 0))
    rng = random.Random(20261006)
    tested = cutoff_nonunit = unit_parity = singular = 0
    for _ in range(3000):
        corrections = tuple(rng.choice(correction_pool) for _ in range(4))
        ok, _, _ = fixture_check(prime_norm_groups, corrections)
        assert ok
        tested += 1
        cutoff_nonunit += 1

    # The even-parity closed case: O_0,...,O_3 are units. Keep the exact
    # formula whenever its Y and determinant hypotheses hold.
    groups_closed = prime_norm_groups[:4] + ((1, 0),) * 4
    ok, _, _ = fixture_check(groups_closed, ((1, 0),) * 4)
    assert ok
    unit_parity += 1

    # Conversely make U_0 real with K_0=conj(B_0); then M>=Q_* and the
    # nonvanishing conclusion must not be asserted.
    raw = sign_products(prime_norm_groups, ((1, 0),) * 4)
    corrections = (gconj(gaussian_product(tuple(
        g if s == 1 else gconj(g) for g, s in zip(prime_norm_groups, SIGN_ROWS[0])))),
                    (1, 0), (1, 0), (1, 0))
    assert gnorm(corrections[0]) == prod(map(gnorm, prime_norm_groups))
    points = sign_products(prime_norm_groups, corrections)
    assert points[0][1] == 0
    try:
        fixture_check(prime_norm_groups, corrections)
    except AssertionError:
        # The strict-cutoff wrapper correctly rejects this deliberately
        # singular fixture. Check the universal cutoff failure directly.
        Gij = pair_g_norms(prime_norm_groups)
        assert max(gnorm(k) for k in corrections) * min(Gij.values()) >= prod(map(gnorm, prime_norm_groups))
        singular += 1
    else:
        raise AssertionError("the constructed real witness unexpectedly passed cutoff")
    return tested, cutoff_nonunit, unit_parity, singular


def walsh_words(m: int) -> tuple[tuple[int, ...], ...]:
    r = 1 << m
    return tuple(tuple(-1 if bin(a & x).count("1") & 1 else 1 for x in range(r))
                 for a in range(r))


def alpha_beta(words: tuple[tuple[int, ...], ...], indices: tuple[int, int, int, int]) -> tuple[int, int]:
    i0, i1, i2, i3 = indices
    alpha = beta = 0
    for x in range(len(words[0])):
        rel = tuple(words[i][x] * words[i0][x] for i in (i1, i2, i3))
        if rel[0] * rel[1] * rel[2] == 1:
            alpha += 1
        else:
            beta += 1
    assert alpha % 4 == beta % 4 == 0
    return alpha // 4, beta // 4


def check_walsh_correlation_counts() -> tuple[tuple[int, int, int, int], ...]:
    results = []
    for m in (2, 3, 4):
        words = walsh_words(m)
        r = len(words[0])
        counts = {}
        closed = 0
        for indices in combinations(range(r), 4):
            ab = alpha_beta(words, indices)
            counts[ab] = counts.get(ab, 0) + 1
            xor = indices[0] ^ indices[1] ^ indices[2] ^ indices[3]
            if xor == 0:
                closed += 1
                assert ab == (r // 4, 0)
            else:
                assert ab == (r // 8, r // 8) if r >= 8 else ab == (0, 0)
        expected_closed = (r * (r - 1) * (r - 2)) // 24
        assert closed == expected_closed
        results.append((r, closed, counts.get((r // 4, 0), 0),
                        counts.get((r // 8, r // 8), 0)))

    # All arithmetic-admissible parity counts for r=4,8,16.
    for r in (4, 8, 16):
        allowed = []
        for alpha in range(r // 4 + 1):
            beta = r // 4 - alpha
            if abs(alpha - beta) * 3 <= alpha + beta:
                allowed.append((alpha, beta))
        assert allowed == ([] if r == 4 else [(r // 8, r // 8)])
    return tuple(results)


def check_independent_walsh_family() -> tuple[int, ...]:
    sizes = []
    for m in range(3, 9):
        r = 1 << m
        coordinate_labels = tuple(1 << i for i in range(m))
        words = walsh_words(m)
        chosen = tuple(words[a] for a in coordinate_labels)
        assert all(sum(x * y for x, y in zip(u, v)) == 0
                   for u, v in combinations(chosen, 2))
        for mask in range(1, 1 << m):
            xor = 0
            for i, label in enumerate(coordinate_labels):
                if (mask >> i) & 1:
                    xor ^= label
            assert xor != 0
        if m >= 3:
            for q in combinations(chosen, 4):
                a, b = alpha_beta(words, tuple(coordinate_labels[chosen.index(v)] for v in q))
                assert a == b == r // 8
        sizes.append(m)
    return tuple(sizes)


def check_exponent_budget() -> tuple[tuple[int, int], ...]:
    outputs = []
    for m in range(2, 6):
        budget = 1 << (m - 1)
        valid = []
        exceptional = {(0,) * m}
        for i in range(m):
            for sign in (-1, 1):
                eps = [0] * m
                eps[i] = sign
                exceptional.add(tuple(eps))
            for j in range(m):
                if i != j:
                    eps = [0] * m
                    eps[i], eps[j] = 1, -1
                    exceptional.add(tuple(eps))
        for eps in product(range(-3, 4), repeat=m):
            cost = sum(abs(sum(eps[i] for i in range(m) if (mask >> i) & 1))
                       for mask in range(1, 1 << m))
            # Exact finite check of the uniform gap in the note: outside zero,
            # signed units, and directed pair differences, F >= 3r/2.
            if eps not in exceptional:
                assert cost >= (3 * budget + 1) // 2
            if cost <= budget:
                valid.append((eps, cost))
        expected = {((0,) * m, 0)}
        for i in range(m):
            for sign in (-1, 1):
                eps = [0] * m
                eps[i] = sign
                expected.add((tuple(eps), budget))
            for j in range(m):
                if i == j:
                    continue
                eps = [0] * m
                eps[i], eps[j] = 1, -1
                expected.add((tuple(eps), budget))
        assert set(valid) == expected
        outputs.append((m, len(valid)))
    return tuple(outputs)


if __name__ == "__main__":
    gaussian_counts = check_gaussian_fixtures()
    code_counts = check_walsh_correlation_counts()
    independent_sizes = check_independent_walsh_family()
    exponent_counts = check_exponent_budget()
    print(f"PASS: {gaussian_counts[0]} generic 8-group Gaussian fixtures; common-factor table, "
          "integer triangle rows, primitive rank-two cross products, and all three finite bounds checked.")
    print(f"PASS: {gaussian_counts[2]} all-even closed fixtures with unit odd groups; "
          f"{gaussian_counts[3]} deliberately singular fixtures fail the strict cutoff as required.")
    print("Walsh four-word counts (r, closed, alpha=r/4 beta=0, balanced): " +
          ", ".join(map(str, code_counts)))
    print(f"PASS: independent Walsh families have {independent_sizes} words for r=8,...,256; "
          "all nonempty products remain nonconstant and all tested quartets balance.")
    print(f"PASS: exact exponent-budget classifications for m=2,...,5: {exponent_counts}.")
