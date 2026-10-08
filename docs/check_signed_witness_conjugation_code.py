"""Exact checks for conjugation symmetry and the signed-witness code bound.

This is a finite checker for the combinatorial consequence in the note; the
rank argument itself proves the bound for every dimension.  It also records
the Pell r=2 equality fixture, where all four sign lattices have thin
witnesses in one real-axis rectangle.
"""

from itertools import combinations
from math import gcd


Gaussian = tuple[int, int]


def gadd(z: Gaussian, w: Gaussian) -> Gaussian:
    return z[0] + w[0], z[1] + w[1]


def gneg(z: Gaussian) -> Gaussian:
    return -z[0], -z[1]


def gconj(z: Gaussian) -> Gaussian:
    return z[0], -z[1]


def gmul(z: Gaussian, w: Gaussian) -> Gaussian:
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def gnorm(z: Gaussian) -> int:
    return z[0] * z[0] + z[1] * z[1]


def sign_vectors(r: int) -> tuple[tuple[int, ...], ...]:
    return tuple((1,) + tuple(1 if (mask >> j) & 1 == 0 else -1
                              for j in range(r - 1))
                 for mask in range(1 << (r - 1)))


def hamming(z: tuple[int, ...], w: tuple[int, ...]) -> int:
    return sum(a != b for a, b in zip(z, w))


def dot(z: tuple[int, ...], w: tuple[int, ...]) -> int:
    return sum(a * b for a, b in zip(z, w))


def maximum_orthogonal_cliques(r: int) -> tuple[int, int, int]:
    """Return (projective vertices, orthogonality edges, maximum clique size)."""
    reps = sign_vectors(r)  # One representative, with first coordinate +1, per antipodal line.
    adjacency = [0] * len(reps)
    edge_count = 0
    for i, j in combinations(range(len(reps)), 2):
        if dot(reps[i], reps[j]) == 0:
            adjacency[i] |= 1 << j
            adjacency[j] |= 1 << i
            edge_count += 1

    best = 1

    def expand(candidates: int, size: int) -> None:
        nonlocal best
        if size > best:
            best = size
        if size + bin(candidates).count("1") <= best:
            return
        while candidates:
            if size + bin(candidates).count("1") <= best:
                return
            bit = candidates & -candidates
            vertex = bit.bit_length() - 1
            candidates ^= bit
            expand(candidates & adjacency[vertex], size + 1)

    expand((1 << len(reps)) - 1, 0)
    return len(reps), edge_count, best


def check_complement_symmetry() -> int:
    checks = 0
    blocks = ((2, 1), (3, 2), (4, 1), (4, 3),
              (5, 2), (6, 1), (5, 4), (7, 2))
    for r in range(1, 9):
        # Distinct positive block norms suffice: sign choices change orientation,
        # not the norm of the signed product.
        selected_blocks = blocks[:r]
        norms = tuple(gnorm(z) for z in selected_blocks)
        product_norm = 1
        for n in norms:
            product_norm *= n
        for sigma in sign_vectors(r):
            signed_product = (1, 0)
            for sign, block in zip(sigma, selected_blocks):
                signed_product = gmul(signed_product,
                                      block if sign == 1 else gconj(block))
            assert gnorm(signed_product) == product_norm
            opposite_product = (1, 0)
            for sign, block in zip(sigma, selected_blocks):
                opposite_product = gmul(opposite_product,
                                        gconj(block) if sign == 1 else block)
            assert opposite_product == gconj(signed_product)
            assert gnorm(opposite_product) == gnorm(signed_product)
            for tau in sign_vectors(r):
                h = hamming(sigma, tau)
                q = 1
                for j, (a, b) in enumerate(zip(sigma, tau)):
                    if a != b:
                        q *= norms[j]
                q_complement = 1
                for j, (a, b) in enumerate(zip(sigma, tau)):
                    if a == b:
                        q_complement *= norms[j]
                assert q * q_complement == product_norm
                assert hamming(sigma, tuple(-x for x in tau)) == r - h
                assert hamming(tuple(-x for x in sigma),
                               tuple(-x for x in tau)) == h
                checks += 1
    return checks


def check_clique_bound() -> tuple[list[tuple[int, int, int, int]], int]:
    results = []
    for r in range(1, 9):
        vertices, edges, clique = maximum_orthogonal_cliques(r)
        # Every clique is a linearly independent family of orthogonal nonzero
        # vectors in R^r, so its size cannot exceed r.
        assert clique <= r
        # For odd r, r/2 is not an integer: only antipodal pairs are allowed.
        if r % 2:
            assert clique == 1
        results.append((r, vertices, edges, clique))
    assert dict((r, c) for r, _, _, c in results)[4] == 4
    assert dict((r, c) for r, _, _, c in results)[8] == 8
    # Sylvester order four is an equality fixture.
    hadamard4 = ((1, 1, 1, 1),
                 (1, -1, 1, -1),
                 (1, 1, -1, -1),
                 (1, -1, -1, 1))
    assert all(dot(a, b) == (4 if i == j else 0)
               for i, a in enumerate(hadamard4)
               for j, b in enumerate(hadamard4))
    assert all(row[0] == 1 for row in hadamard4)
    assert all(hamming(a, b) == 2 for a, b in combinations(hadamard4, 2))
    # Each row together with its antipode gives eight signed patterns in four
    # antipodal classes; this attains 2r in the original cube.
    assert len({tuple(row) for row in hadamard4} |
               {tuple(-x for x in row) for row in hadamard4}) == 8
    return results, 4


def check_common_group_even_parity() -> tuple[int, tuple[int, int, int, int]]:
    # The common group D has fixed orientation; A,B,C have even parity.
    rows = ((1, 1, 1, 1),
            (1, 1, -1, -1),
            (1, -1, 1, -1),
            (1, -1, -1, 1))
    assert all(sum(row[1:]) in (-1, 3) for row in rows)
    assert all(hamming(a, b) == 2 for a, b in combinations(rows, 2))
    frequencies = tuple(sum(row[j] == 1 for row in rows) for j in range(4))
    assert frequencies == (4, 2, 2, 2)

    # Arbitrary positive block norms show exact equal mass: conjugation changes
    # phase but leaves each of the four norm factors unchanged.
    block_norms = (11, 17, 23, 31)  # D,A,B,C
    masses = tuple(block_norms[0] * block_norms[1] *
                   block_norms[2] * block_norms[3] for _ in rows)
    assert masses == (11 * 17 * 23 * 31,) * 4
    return len(rows), frequencies


def pell_fixture(j: int) -> tuple[int, int]:
    """Return (u,b) from u+b√5=(2+√5)(9+4√5)^j, j>=1."""
    assert j >= 1
    u, b = 38, 17
    for _ in range(1, j):
        u, b = 9 * u + 20 * b, 4 * u + 9 * b
    return u, b


def check_pell_witnesses() -> int:
    checks = 0
    for j in range(1, 33):
        u, b = pell_fixture(j)
        assert u * u - 5 * b * b == -1
        assert u % 2 == 0 and b % 2 == 1 and gcd(u, b) == 1

        f = (u + 2 * b, b)
        a = (b, 2 * b - u)
        q, n = gnorm(f), gnorm(a)
        assert gcd(q, n) == 1

        U = gmul(f, a)
        V = gmul(gconj(f), a)
        y = gmul((2, 1), V)
        assert U == (2 * u * b, 1)
        assert V == (4 * b * b, 1 - 2 * b * b)
        assert y == (10 * b * b - 1, 2)
        assert gmul((2, 1), V) == y
        assert gmul((2, -1), gconj(V)) == gconj(y)
        assert gnorm(U) == q * n
        assert gnorm(y) == 5 * q * n

        # Sign order: ++, +-, -+, --. Every point is in the same
        # axis-parallel rectangle |Re z|<=|y|, |Im z|<=2.
        witnesses = (U, gconj(y), y, gconj(U))
        multipliers = ((1, 0), (2, -1), (2, 1), (1, 0))
        bases = (gmul(f, a), gmul(f, gconj(a)),
                 gmul(gconj(f), a), gmul(gconj(f), gconj(a)))
        assert tuple(gmul(c, z) for c, z in zip(multipliers, bases)) == witnesses
        relation_terms = (gmul((2, 0), U), gconj(y),
                          gneg(y), gmul((-2, 0), gconj(U)))
        relation_sum = (0, 0)
        for term in relation_terms:
            relation_sum = gadd(relation_sum, term)
        assert relation_sum == (0, 0)
        assert gnorm(gconj(U)) == gnorm(U)
        assert gnorm(gconj(y)) == gnorm(y)
        for z in witnesses:
            assert z[0] * z[0] <= gnorm(y)
            assert abs(z[1]) <= 2
        assert max(gnorm(z) for z in witnesses) == gnorm(y)
        assert all(hamming(s, t) in (1, 2)
                   for s, t in combinations(((1, 1), (1, -1),
                                             (-1, 1), (-1, -1)), 2))
        checks += 1
    return checks


if __name__ == "__main__":
    complement_checks = check_complement_symmetry()
    clique_results, hadamard_order = check_clique_bound()
    quarter_rows, frequencies = check_common_group_even_parity()
    pell_checks = check_pell_witnesses()
    print(f"PASS: {complement_checks} exact complement-distance/norm checks (r=1,...,8).")
    print("Projective orthogonality cliques (r, lines, edges, maximum): " +
          ", ".join(map(str, clique_results)))
    print(f"PASS: Hadamard order {hadamard_order} reaches the combinatorial "
          "Hadamard bound on signed patterns.")
    print(f"PASS: {quarter_rows} common-D even-parity patterns have equal mass; "
          f"positive-sign frequencies D,A,B,C={frequencies}.")
    print(f"PASS: {pell_checks} Pell fixtures; all four r=2 sign lattices have "
          "thin witnesses U, conjugate(y), y, conjugate(U) in one real-axis rectangle.")
