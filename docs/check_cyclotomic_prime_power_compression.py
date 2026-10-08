"""Exact finite checks of the prime-power compression in cyclotomic_uniform_rank.md.

Ranks are computed over Q by integer row operations, with row-content
division only after an exact elimination step. No finite-field or
floating-point rank is used. All supports in the 3^a 5^b box, 0<=a,b<=2,
are checked, together with deterministic three-prime and deeper-power
fixtures. These finite checks supplement the general proof.
All support and selected-row subsets of a six-element divisor downset
also verify the arbitrary-row extension, including non-interval row sets.
"""

from __future__ import annotations

from functools import cache
from math import gcd


@cache
def factors(n: int) -> tuple[tuple[int, int], ...]:
    result = []
    p = 2
    while p * p <= n:
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        if exponent:
            result.append((p, exponent))
        p += 1
    if n > 1:
        result.append((n, 1))
    return tuple(result)


@cache
def divisors(n: int) -> tuple[int, ...]:
    result = [1]
    for p, exponent in factors(n):
        result = [d * p**j for d in result for j in range(exponent + 1)]
    return tuple(sorted(result))


@cache
def phi(n: int) -> int:
    answer = n
    for p, _ in factors(n):
        answer = answer // p * (p - 1)
    return answer


@cache
def mu(n: int) -> int:
    fs = factors(n)
    return 0 if any(a > 1 for _, a in fs) else (-1) ** len(fs)


def charged_cost(support: frozenset[int]) -> int:
    return sum(2 if e == 1 else phi(e) for e in support)


def rational_rank(matrix: list[list[int]], columns: int) -> int:
    """Fraction-free elimination; the allowed row operations preserve Q-rank."""
    matrix = [row[:] for row in matrix]
    rank = 0
    for col in range(columns):
        pivot = next((j for j in range(rank, len(matrix)) if matrix[j][col]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        value = matrix[rank][col]
        for j in range(rank + 1, len(matrix)):
            coefficient = matrix[j][col]
            if not coefficient:
                continue
            row = [value * x - coefficient * y for x, y in zip(matrix[j], matrix[rank])]
            content = 0
            for x in row:
                content = gcd(content, x)
            matrix[j] = [x // content for x in row] if content else row
        rank += 1
        if rank == len(matrix):
            break
    return rank


@cache
def nullity(support: frozenset[int], h: int) -> int:
    columns = sorted(support)
    # Every other low row is identically zero because it divides no column.
    rows = sorted({d for e in columns for d in divisors(e) if d < h})
    matrix = [[mu(e // d) if e % d == 0 else 0 for e in columns] for d in rows]
    return len(columns) - rational_rank(matrix, len(columns))


def compression_stages(support: frozenset[int], p: int):
    slices: dict[int, set[int]] = {}
    for e in support:
        a, s = 0, e
        while s % p == 0:
            a += 1
            s //= p
        slices.setdefault(a, set()).add(s)
    for a in range(max(slices, default=0), 0, -1):
        top = set(slices.get(a, set()))
        common = set(top)
        for j in range(a):
            common.intersection_update(slices.get(j, set()))
        assert slices.get(a + 1, set()) <= common
        for j in range(a):
            slices.setdefault(j, set()).update(top)
        slices[a] = common
        yield frozenset(p**j * s for j, values in slices.items() for s in values)


def check_support(support: frozenset[int], thresholds: tuple[int, ...]) -> int:
    original = support
    comparisons = 0
    for p in sorted({p for e in support for p, _ in factors(e)}):
        for compressed in compression_stages(support, p):
            assert charged_cost(compressed) <= charged_cost(support)
            for h in thresholds:
                assert nullity(compressed, h) >= nullity(support, h), (support, p, h)
                comparisons += 1
            support = compressed
    assert all(set(divisors(e)) <= support for e in support)
    assert charged_cost(support) <= charged_cost(original)
    for h in thresholds:
        assert nullity(support, h) == sum(e >= h for e in support)
    return comparisons


@cache
def selected_nullity(support: frozenset[int], rows: frozenset[int]) -> int:
    columns = sorted(support)
    matrix = [[mu(e // d) if e % d == 0 else 0 for e in columns]
              for d in sorted(rows)]
    return len(columns) - rational_rank(matrix, len(columns))


def check_arbitrary_selected_rows() -> tuple[int, int]:
    universe = (1, 3, 5, 9, 15, 45)
    subsets = [frozenset(e for j, e in enumerate(universe) if mask & (1 << j))
               for mask in range(1 << len(universe))]
    cases = comparisons = 0
    for original in subsets:
        for rows in subsets:
            support = original
            for p in sorted({p for e in support for p, _ in factors(e)}):
                for compressed in compression_stages(support, p):
                    assert charged_cost(compressed) <= charged_cost(support)
                    assert sum(map(phi, compressed)) <= sum(map(phi, support))
                    assert selected_nullity(compressed, rows) >= selected_nullity(support, rows), (
                        support, compressed, rows, p
                    )
                    comparisons += 1
                    support = compressed
            assert all(set(divisors(e)) <= support for e in support)
            assert selected_nullity(support, rows) == len(support - rows)
            cases += 1
    return cases, comparisons


def main() -> None:
    # Ordinary single-column shifting would fail on this example.
    assert nullity(frozenset({9}), 3) == 1
    assert nullity(frozenset({3}), 3) == 0
    assert list(compression_stages(frozenset({9}), 3))[-1] == frozenset({1, 3})
    assert nullity(frozenset({1, 3}), 3) == 1
    assert charged_cost(frozenset({1, 3})) < charged_cost(frozenset({9}))

    thresholds = tuple(range(1, 32, 2))
    universe = tuple(3**a * 5**b for a in range(3) for b in range(3))
    total = 0
    for mask in range(1 << len(universe)):
        support = frozenset(e for j, e in enumerate(universe) if mask & (1 << j))
        total += check_support(support, thresholds)

    three_prime = tuple(3**a * 5**b * 7**c for a in range(3) for b in range(2) for c in range(2))
    for seed in range(1, 129):
        mask = (seed * 1093) % (1 << len(three_prime))
        support = frozenset(e for j, e in enumerate(three_prime) if mask & (1 << j))
        total += check_support(support, (1, 3, 5, 11, 31, 101, 501))

    for support in (
        frozenset({3**7}),
        frozenset({1, 3**7}),
        frozenset({3**5, 3**3 * 5**2, 3 * 5**4, 5**5}),
        frozenset({3**4 * 5**3, 3**2 * 7**4, 5**2 * 7**3}),
    ):
        total += check_support(support, (1, 3, 9, 27, 101, 501, 1501))

    print(
        f"PASS: {total} exact rational-nullity comparisons; charged cost never "
        "increases, and every final support is divisor closed."
    )
    cases, comparisons = check_arbitrary_selected_rows()
    print(
        f"PASS: {cases} arbitrary support/row-set pairs; {comparisons} exact "
        "compression comparisons; final nullity equals the omitted-column count."
    )


if __name__ == "__main__":
    main()
