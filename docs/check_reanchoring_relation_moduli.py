"""Exact finite checks for reanchoring_relation_moduli.md (stdlib only)."""
from itertools import permutations, product


def check_support():
    m = 4
    errors = [(0, 0, 0, 0), (1, 0, -1, 0), (-1, 1, 0, 1)]
    count = 0
    for old in product((0, 1), repeat=m):
        for new in product((0, 1), repeat=m):
            if len(set(new)) == 1:
                continue
            compatible = new == old or new == tuple(1 - x for x in old)
            for e in range(4):
                for f in range(1, 5):
                    for delta in errors:
                        for c in range(-3, 4):
                            for sign, orientation in product((-1, 1), repeat=2):
                                actual = [c + sign * (e * old[i] + delta[i]) for i in range(m)]
                                remainder = [actual[i] - orientation * f * new[i] for i in range(m)]
                                cost = sum(map(abs, delta)) + sum(map(abs, remainder))
                                assert f <= (e if compatible else 0) + cost
                                count += 1
    print(f"Signed prime-power support inequality: {count} exact cases pass.")


def action(bits, perm):
    full = (0,) + tuple(bits)
    return tuple(full[perm[j]] ^ full[perm[0]] for j in range(1, len(full)))


def check_boolean_group():
    for m in range(2, 5):
        vectors = list(product((0, 1), repeat=m))
        group = {tuple(action(v, p) for v in vectors) for p in permutations(range(m + 1))}
        factorial = 1
        for x in range(1, m + 2):
            factorial *= x
        assert len(group) == factorial
        ones = (1,) * m
        for a in range(1, m + 1):
            perm = list(range(m + 1))
            perm[0], perm[a] = perm[a], perm[0]
            assert action(ones, perm) == tuple(int(j == a) for j in range(1, m + 1))
    orbit = {action((1, 1, 0, 0, 0, 0), p) for p in permutations(range(7))}
    assert len(orbit) == 21
    assert {sum(b) for b in orbit} == {2, 5}
    print("Faithful anchor-permutation group and six-row 21-cut orbit pass.")


def bracket(x, y):
    return x[0] * y[1] - x[1] * y[0]


def invariant(rows):
    a, b, c, d = rows
    return (4 * bracket(a, b) * bracket(c, d) - bracket(a, c) * bracket(b, d)) ** 2


def check_replacement():
    rows = [(2 * j, 1) for j in range(1, 5)]
    assert invariant(rows) == 0
    rows[0] = (1, 0)
    assert invariant(rows) == 16
    print("Replacing an active row by the old anchor destroys the exact invariant zero.")


if __name__ == "__main__":
    check_support()
    check_boolean_group()
    check_replacement()
