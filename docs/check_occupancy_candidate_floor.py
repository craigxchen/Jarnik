"""Exact checks for the exploratory occupancy count M(n,d).

The checker establishes the combinatorics and finite comparisons only. It
does not prove that M is an arithmetic coefficient floor.
"""

from math import factorial

from check_higher_rank_terminal_fusion_grid import case


def occupancy_formula(n, d):
    q = 1 << d
    total = 0
    for full in range(n + 1):
        for empty in range(n - full + 1):
            partial = n - full - empty
            ways = factorial(n) // (
                factorial(full) * factorial(empty) * factorial(partial)
            )
            total += ways * (q - 2) ** partial * min(full, empty)
    return total


def occupancy_by_subsets(n, d):
    block_size = d
    m = n * d
    total = 0
    for subset in range(1 << m):
        full = empty = 0
        for block in range(n):
            block_mask = ((1 << block_size) - 1) << (block * block_size)
            selected = subset & block_mask
            if selected == block_mask:
                full += 1
            elif selected == 0:
                empty += 1
        total += min(full, empty)
    return total


def main():
    for n, d in ((2, 2), (3, 2), (4, 2), (5, 2), (2, 3), (3, 3), (4, 3)):
        assert occupancy_formula(n, d) == occupancy_by_subsets(n, d)

    grid = {}
    for n in range(2, 9):
        for d in range(2, 13):
            item = case(n, d)
            candidate = occupancy_formula(n, d)
            grid[n, d] = item, candidate
            assert item["B"] >= item["h"] * candidate
            assert (item["B"] == item["h"] * candidate) == (d == 2)

    for d, n_max in ((3, 12), (4, 12), (5, 12), (6, 12)):
        for n in range(2, n_max + 1):
            item = case(n, d)
            assert item["B"] > item["h"] * occupancy_formula(n, d)

    print("Subset enumeration agrees with the multinomial formula.")
    print("77 grid cases: equality exactly for d=2; otherwise B/h > M.")
    print("Additional d=3..6 samples through n=12: B/h > M throughout.")


if __name__ == "__main__":
    main()
