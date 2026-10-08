"""Exact audit of the four-row product reduction; no moment solution claimed."""

from fractions import Fraction
from itertools import combinations


def mul(p, q):
    out = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return out


def product(factors):
    out = [Fraction(1)]
    for factor in factors:
        out = mul(out, factor)
    return out


def scale(p, c):
    return [c * x for x in p]


def sub(p, q):
    return [
        (p[i] if i < len(p) else 0) - (q[i] if i < len(q) else 0)
        for i in range(max(len(p), len(q)))
    ]


def audit_four_row_factors():
    # Walsh order eight, minus its constant and one balanced column.
    full = [
        [1 if bin(i & j).count("1") % 2 == 0 else -1 for j in range(8)]
        for i in range(8)
    ]
    b = [row[1] for row in full]
    S = [[row[j] for j in range(2, 8)] for row in full]
    a = list(range(1, 7))
    Arows = [i for i in range(8) if b[i] == 1]
    Brows = [i for i in range(8) if b[i] == -1]

    def G(i, j):
        return product(
            [[Fraction(-S[i][h] * a[h]), Fraction(1)]
             for h in range(6) if S[i][h] == S[j][h]]
        )

    checks = 0
    for i, j in combinations(Arows, 2):
        for k, l in combinations(Brows, 2):
            classes = {"A": [], "B": [], "C": []}
            common = []
            for h in range(6):
                si, sj, sk, sl = (S[r][h] for r in (i, j, k, l))
                if si == sj == sk == sl:
                    common.extend([[Fraction(-si * a[h]), Fraction(1)]] * 2)
                elif si == sk and sj == sl and si == -sj:
                    classes["A"].append(h)
                elif si == sl and sj == sk and si == -sj:
                    classes["B"].append(h)
                elif si == sj and sk == sl and si == -sk:
                    classes["C"].append(h)
                else:
                    majority = 1 if si + sj + sk + sl > 0 else -1
                    common.append([Fraction(-majority * a[h]), Fraction(1)])

            assert len(classes["A"]) == len(classes["B"])
            assert len(classes["A"]) == len(classes["C"]) + 1
            H = product(common)

            def F(name):
                return product(
                    [[Fraction(-a[h] * a[h]), Fraction(0), Fraction(1)]
                     for h in classes[name]]
                )

            assert mul(G(i, k), G(j, l)) == mul(H, F("A"))
            assert mul(G(i, l), G(j, k)) == mul(H, F("B"))
            assert mul([Fraction(0), Fraction(0), Fraction(1)],
                       mul(G(i, j), G(k, l))) == mul(
                           H, mul([Fraction(0), Fraction(0), Fraction(1)], F("C"))
                       )
            checks += 1
    assert checks == 36
    return checks


def audit_local_feasibility():
    # A positive-node solution of the reduced four-row identity only.
    # These five nodes are not asserted to come from an eight-row solution.
    FA = product([[-1, 1], [-4, 1]])
    FB = product([[-2, 1], [-Fraction(8, 3), 1]])
    FC = [-6, 1]
    assert sub(scale(FA, 4), scale(FB, 3)) == mul([0, 1], FC)
    assert Fraction(4, 3) == Fraction(2 * Fraction(8, 3), 1 * 4)
    QA = product([[1, 1], [1, Fraction(1, 4)]])
    QB = product([[1, Fraction(1, 2)], [1, Fraction(3, 8)]])
    d1, d2 = sub(QA, QB)[1:]
    assert (d1, d2) == (Fraction(3, 8), Fraction(1, 16))
    assert d1 / d2 == 6
    assert order_parity_holds([1, 4], [2, 3], [5], 1)
    assert not order_parity_holds([1, 3], [2, 4], [5], 1)
    assert not order_parity_holds([1, 3], [2, 4], [5], -1)


def order_parity_holds(A, B, C, epsilon):
    def above(E, h):
        return sum(e > h for e in E)

    return (
        all(epsilon == (-1) ** (1 + above(B, h) + above(C, h)) for h in A)
        and all(epsilon == (-1) ** (above(A, h) + above(C, h)) for h in B)
        and all((above(A, h) + above(B, h)) % 2 == 0 for h in C)
    )


def audit_admissible_prefix_example():
    # The published 11-column terminal path for the even-s, q=2 automaton.
    # It lacks four required distance increments, so is not a full template.
    word = [128, 64, 206, 205, 35, 19, 240, 8, 4, 60, 195]
    i, j, k, l = 0, 2, 4, 6
    A, B, C = [], [], []
    for h, mask in enumerate(word, 1):
        si, sj, sk, sl = (1 if mask & (1 << r) else -1 for r in (i, j, k, l))
        if si == sk and sj == sl and si == -sj:
            A.append(h)
        elif si == sl and sj == sk and si == -sj:
            B.append(h)
        elif si == sj and sk == sl and si == -sk:
            C.append(h)
    assert (A, B, C) == ([3, 6], [10, 11], [7])
    assert all(a < b for a, b in zip(A, B))
    assert order_parity_holds(A, B, C, 1)
    rank = {row: pos for pos, row in enumerate((6, 7, 0, 1, 2, 3, 4, 5))}
    v = (rank[i] - rank[l]) * (rank[j] - rank[k])
    w = (rank[i] - rank[j]) * (rank[k] - rank[l])
    assert v * w > 0


def audit_two_quartet_geometry():
    # A rational spectral ordering fixture for the exact-value bound.
    c = {6: 0, 7: 1, 0: 2, 1: 3, 2: 4, 3: 5, 4: 6, 5: 7}

    def ratio(i, j, k, l):
        u = (c[i] - c[k]) * (c[j] - c[l])
        v = (c[i] - c[l]) * (c[j] - c[k])
        return Fraction(u, v)

    low = ratio(1, 2, 6, 7)
    high = ratio(0, 1, 4, 5)
    assert low == Fraction(9, 8)
    assert high == Fraction(16, 15)
    assert (low - 1) * (high - 1) == Fraction(1, 120) < 1


if __name__ == "__main__":
    count = audit_four_row_factors()
    audit_local_feasibility()
    audit_admissible_prefix_example()
    audit_two_quartet_geometry()
    print(f"checked {count} Walsh cross quadruples, rational identities, and an admissible prefix example")
