#!/usr/bin/env python3
"""Exact finite checks for the odd-moment degree-window note.

This is a standard-library checker: all arithmetic is integral and the
small Newton scan is deliberately bounded.  It checks the collision
polynomial, the defect-window dynamic program, the padded incidence
identities, and the two Hadamard fixtures used as sanity checks.
"""

from itertools import combinations, product
from fractions import Fraction


def powers(xs, k):
    return sum(x**k for x in xs)


def elementary_from_powers(xs):
    """Newton recurrence, retained as Fractions to make exactness explicit."""
    d = len(xs)
    p = [Fraction(0)] + [Fraction(powers(xs, k)) for k in range(1, d + 1)]
    e = [Fraction(0)] * (d + 1)
    e[0] = Fraction(1)
    for k in range(1, d + 1):
        e[k] = sum((-1) ** (i - 1) * e[k - i] * p[i] for i in range(1, k + 1)) / k
    return e


def check_newton_collision():
    # A bounded exact sanity scan.  The general argument is the Newton
    # identity in the accompanying note; this checks every nonzero tuple
    # in {-3,...,-1,1,...,3} for d <= 6.
    vals = (-3, -2, -1, 1, 2, 3)
    for d in range(1, 7):
        for xs in product(vals, repeat=d):
            if any(powers(xs, k) for k in (1, 3, 5) if k <= d):
                continue
            e = elementary_from_powers(xs)
            assert all(e[k] == 0 for k in range(1, d + 1, 2))
            assert d % 2 == 0
            assert tuple(sorted(xs)) == tuple(sorted(-x for x in xs))
    # The strict support threshold matters: this five-term relation
    # cancels orders one and three without an opposite pair.
    sharp = (1, 5, 9, -7, -8)
    assert powers(sharp, 1) == powers(sharp, 3) == 0
    assert len(set(map(abs, sharp))) == 5 and powers(sharp, 5) != 0


def scalar_max_excess(m):
    s = (m + 1) // 4
    budget = 16 * m - 56 * s - 40
    if budget < 0:
        return None
    best = None
    for n1 in range(m + 1):
        for n2 in range(m - n1 + 1):
            for n3 in range(m - n1 - n2 + 1):
                n4 = m - n1 - n2 - n3
                defect = 9 * n1 + 4 * n2 + n3
                if defect <= budget:
                    excess = defect + n4 - 2 * m
                    best = excess if best is None else max(best, excess)
    return best


def check_defect_windows():
    assert [scalar_max_excess(m) for m in (6, 10, 13, 14, 17, 18)] == [
        -6,
        -4,
        -13,
        -2,
        -11,
        1,
    ]
    assert 16 * 19 - 56 * 5 - 40 < 0
    assert scalar_max_excess(20) == -20
    assert scalar_max_excess(21) == -9


def check_degree_eighteen_profiles():
    positive, zero = [], []
    for n1 in range(19):
        for n2 in range(19 - n1):
            for n3 in range(19 - n1 - n2):
                n4 = 18 - n1 - n2 - n3
                defect = 9 * n1 + 4 * n2 + n3
                excess = defect + n4 - 36
                if excess < 0:
                    continue
                for beta in range(5):
                    D = defect + beta * beta + 16
                    B = 3 * n1 + 2 * n2 + n3 + beta + 4
                    K, q = Fraction(40 - D, 2), Fraction(2 * B - 20, 4)
                    if K < 0 or K.denominator != 1 or q.denominator != 1:
                        continue
                    fixture = (n1, n2, n3, n4, beta, K, q)
                    if excess > 0:
                        assert n1 and q + K < 3
                        positive.append(fixture)
                    else:
                        assert n1 == n3 == beta == K == 0
                        assert n2 == 6 and n4 == 12 and q == 3
                        # Every W_i is even, but K=0 forces W_i=q=3.
                        zero.append(fixture)
    assert positive == [
        (2, 1, 0, 15, 0, 1, 1),
        (2, 1, 1, 14, 1, 0, 2),
        (2, 1, 2, 13, 0, 0, 2),
    ]
    assert zero == [(0, 6, 0, 12, 0, 0, 3)]
    assert [m for m in range(1, 22) if scalar_max_excess(m) is not None] == [
        6, 10, 13, 14, 17, 18, 20, 21,
    ]
    assert all(scalar_max_excess(m) < 0
               for m in (6, 10, 13, 14, 17, 20, 21))


def orient_columns(rows):
    columns = []
    for j in range(len(rows[0])):
        col = tuple(row[j] for row in rows)
        if sum(x == 1 for x in col) < 4:
            col = tuple(-x for x in col)
        minority = frozenset(i for i, x in enumerate(col) if x == -1)
        b = 4 - len(minority)
        columns.append((col, minority, b))
    return columns


def check_padded_identities(rows):
    columns = orient_columns(rows)
    L = len(columns)
    D = sum(b * b for _, _, b in columns)
    B = sum(b for _, _, b in columns)
    distances = {}
    for i, j in combinations(range(8), 2):
        d = sum(col[i] != col[j] for col, _, _ in columns)
        distances[i, j] = d
        assert d % 2 == 0 and d >= L // 2
    A = {}
    for i, j in combinations(range(8), 2):
        value = (2 * distances[i, j] - L)
        assert value % 4 == 0 and value >= 0
        A[i, j] = A[j, i] = value // 4
    K = sum(A[i, j] for i, j in combinations(range(8), 2))
    assert K == L - Fraction(D, 2)
    q_num = 2 * B - L
    assert q_num % 4 == 0
    q = q_num // 4
    for i in range(8):
        W = sum(b for _, minority, b in columns if i in minority)
        degree = sum(A.get((i, j), 0) for j in range(8) if j != i)
        assert W == q + degree


def paley_hadamard_12():
    q = 11
    legendre = {
        x: 0 if x % q == 0 else (1 if pow(x % q, 5, q) == 1 else -1)
        for x in range(q)
    }
    core = [
        [(-1 if a == b else legendre[(a - b) % q]) for b in range(q)]
        for a in range(q)
    ]
    return [[1] * 12] + [[1] + row for row in core]


def check_hadamard_fixtures():
    H12 = paley_hadamard_12()
    assert all(sum(H12[i][k] * H12[j][k] for k in range(12)) == (12 if i == j else 0)
               for i in range(12) for j in range(12))

    # An H12 puncture: eight rows, with two columns deleted, has d_min=5.
    punctured = [[H12[i][j] for j in range(12) if j not in (0, 1)] for i in range(8)]
    assert min(sum(a != b for a, b in zip(x, y)) for x, y in combinations(punctured, 2)) == 5
    assert sorted(len(minority) for _, minority, _ in orient_columns(punctured)).count(3) == 7

    H2 = ((1, 1), (1, -1))
    chosen = (0, 1, 2, 4)
    tensor = [
        [H12[i][j] * H2[e][delta] for j in range(12) for delta in range(2)]
        for i in chosen for e in range(2)
    ]
    # Delete the constant and one balanced tensor column.
    tensor22 = [[x for k, x in enumerate(row) if k not in (0, 1)] for row in tensor]
    assert min(sum(a != b for a, b in zip(x, y)) for x, y in combinations(tensor22, 2)) == 11
    profile = [len(minority) for _, minority, _ in orient_columns(tensor22)]
    assert profile.count(2) == 8 and profile.count(4) == 14


def main():
    check_newton_collision()
    check_defect_windows()
    check_degree_eighteen_profiles()
    # H16 rows indexed by a 3-dimensional subspace of F_2^4, plus the
    # eight even-parity words of length four.  This gives L=20 and tests
    # the nonzero-A padded identities using exact signs only.
    baseline = [
        [1 if (bin(r & c).count("1") % 2 == 0) else -1 for c in range(16)]
        for r in range(8)
    ]
    extra = []
    for u in range(8):
        bits = [(u >> k) & 1 for k in range(3)]
        bits.append(sum(bits) % 2)
        extra.append([1 if bit == 0 else -1 for bit in bits])
    check_padded_identities([row + extra[i] for i, row in enumerate(baseline)])
    check_hadamard_fixtures()
    print("odd-moment degree-window checks: OK")


if __name__ == "__main__":
    main()
