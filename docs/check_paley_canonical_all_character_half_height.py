"""Exact finite audit of the canonical Paley all-character theorem."""

from itertools import product
from random import Random

from check_general_four_row_quartic_slack_obstruction import paley


def physical_matrix(h, b):
    m, q = len(h), len(h) - 1
    r = b * q
    s = [[h[x][j % q + 1] for j in range(r)] for x in range(m)]
    for x in range(1, m):
        s[x][x - 1] *= -1  # Finite row x uses its diagonal label x.
    s[0][q] *= -1         # Infinity uses a second copy of label 1.
    return s


def audit_vector(h, s, c, direct):
    m, b = len(h), len(s[0]) // (len(h) - 1)
    r = b * (m - 1)
    t = [sum(c[x] * h[x][a] for x in range(m)) for a in range(m)]
    assert t[0] == 0
    l1 = sum(abs(v) for v in t)
    support = sum(v != 0 for v in c)
    assert l1 >= m
    assert sum(v * v for v in t) == m * support
    beneficial_finite = sum(c[x] * t[x] < 0 for x in range(1, m) if c[x])
    beneficial_infinity = int(c[0] * t[1] > 0)
    k = beneficial_finite + beneficial_infinity
    assert sum(c[x] * t[x] for x in range(1, m)) == -support
    if support >= 4:
        defect = sum(abs(v) * (support - abs(v)) for v in t)
        assert defect == support * (l1 - m)
        # H_(x,x)=-1, so c_x T_x<0 is the beneficial condition.
        perfect_beneficial = sum(
            c[x] * t[x] < 0 and abs(t[x]) == support
            for x in range(1, m) if c[x]
        )
        assert perfect_beneficial <= 1
        assert defect >= 2 * (support - 2) * max(0, k - 2)
    else:
        assert support == 2 and k <= 2
    assert b * l1 % 2 == 0
    budget = b * l1 // 2 + support - 2 * k
    assert budget >= b * m // 2 - 2
    if direct:
        literal = 0
        for j in range(r):
            value = sum(c[x] * s[x][j] for x in range(m))
            assert value % 2 == 0
            literal += abs(value) // 2
        assert literal == budget
    return budget


def audit_order(q, exhaustive):
    h, _ = paley(q)
    m, b = len(h), 5
    assert all(sum(h[x][a] * h[y][a] for a in range(m))
               == (m if x == y else 0)
               for x in range(m) for y in range(m))
    d = [1] + [-1] * q
    k_matrix = [[d[x] * h[x][a] for a in range(m)] for x in range(m)]
    assert all(k_matrix[x][x] == 1 for x in range(m))
    assert all(k_matrix[x][y] == -k_matrix[y][x]
               for x in range(m) for y in range(x + 1, m))
    s = physical_matrix(h, b)
    rng = Random(156738)
    if exhaustive:
        vectors = (c for c in product((-1, 0, 1), repeat=m)
                   if any(c) and sum(c) == 0)
    else:
        def samples():
            for _ in range(20000):
                size = rng.randrange(1, m // 2 + 1)
                rows = rng.sample(range(m), 2 * size)
                c = [0] * m
                for x in rows[:size]:
                    c[x] = 1
                for x in rows[size:]:
                    c[x] = -1
                yield tuple(c)
        vectors = samples()
    count, minimum = 0, None
    for c in vectors:
        budget = audit_vector(h, s, c, direct=(not exhaustive or count % 100 == 0))
        minimum = budget if minimum is None else min(minimum, budget)
        count += 1
    equality = [0] * m
    equality[0], equality[1] = 1, -1
    assert audit_vector(h, s, equality, True) == b * m // 2 - 2
    assert minimum == b * m // 2 - 2
    print(f"M={m}: {count} ternary vectors, exact min B={minimum}")


def main():
    audit_order(11, True)
    audit_order(19, False)
    audit_order(23, False)
    print("PASS: skew-Paley identities and sharp canonical all-character minimum")


if __name__ == "__main__":
    main()
