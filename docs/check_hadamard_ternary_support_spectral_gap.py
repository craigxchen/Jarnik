"""Finite exact audit of the high-support Hadamard transform gap."""

from itertools import product

from check_general_four_row_quartic_slack_obstruction import paley


def walsh(m):
    return [
        [1 if bin(x & a).count("1") % 2 == 0 else -1
         for a in range(m)]
        for x in range(m)
    ]


def audit(h_matrix):
    m, b = len(h_matrix), 5
    checked, min_gap = 0, None
    for c in product((-1, 0, 1), repeat=m):
        if sum(c) != 0 or all(x == 0 for x in c):
            continue
        h = sum(x != 0 for x in c)
        transform = [sum(c[x] * h_matrix[x][a] for x in range(m))
                     for a in range(m)]
        l1 = sum(abs(x) for x in transform)
        s = sum(x != 0 for x in transform)
        delta_m = l1 - m
        assert delta_m >= 0
        assert 2 * h * s <= (h + 2) * l1 - h * m
        full_column = max(abs(x) for x in transform) == m
        if not full_column:
            assert l1 >= 2 * h
            min_gap = l1 - 2 * h if min_gap is None else min(min_gap, l1 - 2 * h)
        else:
            assert h == m
        # Put the flip at row x in physical column x; its base label
        # is 1+(x mod (m-1)), so each label receives at most two flips.
        base_b = b * l1 // 2
        assert b * l1 % 2 == 0
        correction = 0
        beneficial = 0
        for x in range(m):
            a = 1 + x % (m - 1)
            coefficient = transform[a] // 2
            assert transform[a] % 2 == 0
            correction += abs(coefficient - c[x] * h_matrix[x][a]) - abs(coefficient)
            if c[x] and coefficient * c[x] * h_matrix[x][a] > 0:
                beneficial += 1
        physical_b = base_b + correction
        assert physical_b == base_b + h - 2 * beneficial
        assert beneficial <= 2 * s
        # Equation (6), multiplied by 2h to retain integer arithmetic.
        assert (2 * physical_b - b * (m - 1)) * h >= (
            b * h + 2 * h * h - 8 * m
            + delta_m * ((b - 4) * h - 8)
        )
        if not full_column:
            assert physical_b >= (b - 1) * h
            if 2 * (b - 1) * h > b * (m - 1):
                assert 2 * physical_b > b * (m - 1)
        else:
            assert 2 * physical_b > b * (m - 1)
        if m >= 12 and h >= 8 and h * h >= 4 * m:
            assert 2 * physical_b > b * (m - 1)
        checked += 1
    print(f"M={m}: {checked} zero-sum ternary characters, min L1-2h={min_gap}")


def main():
    audit(walsh(8))
    audit(paley(11)[0])
    print("PASS: universal transform gap and five-copy high-support budget")


if __name__ == "__main__":
    main()
