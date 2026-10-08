"""Exhaustive small-budget check for the truncated cyclotomic Möbius matrix.

For each odd h, enumerate every finite support E of odd indices with

    sum(phi(e) for e in E, e > 1) + 2*[1 in E] <= 4*h.

The coefficient 2 at e=1 is the minimum nonzero width compatible with
the even base-layer parity of the endpoint fibre. On odd d<h and e in E,
the matrix entry is mu(e/d) when d|e, and zero otherwise. Its nullity is
computed exactly over F_101. Since rank mod 101 is at most rational rank,
an upper bound for modular nullity is also an upper bound over Q.

The candidate search is complete: phi(e)^2 >= e for every odd e. Thus
phi(e)<=4*h implies e<=(4*h)^2, the sieve limit used below.

Default: odd h<=19, quick check. Use --full for odd h<=31 (about one
minute on the original host), or --max-h N for another odd limit.
The universal bound of three is FALSE. The exact divisor-closed
counterexamples in cyclotomic_nullity_four_counterexample.md have
nullity four within the charged budget. They are checked separately
before the small-threshold enumeration; use --counterexamples-only
to run just those checks. The exhaustive h<=31 results remain valid.
"""

from __future__ import annotations

import argparse
from functools import cache
from itertools import combinations


MODULUS = 101


@cache
def factorization(n: int) -> tuple[tuple[int, int], ...]:
    factors = []
    p = 2
    while p * p <= n:
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        if exponent:
            factors.append((p, exponent))
        p += 1
    if n > 1:
        factors.append((n, 1))
    return tuple(factors)


@cache
def divisors(n: int) -> tuple[int, ...]:
    result = [1]
    for p, exponent in factorization(n):
        result = [d * p**j for d in result for j in range(exponent + 1)]
    return tuple(sorted(result))


def exact_phi(n: int) -> int:
    result = n
    for p, _ in factorization(n):
        result = result // p * (p - 1)
    return result


def exact_mu(n: int) -> int:
    factors = factorization(n)
    return 0 if any(a > 1 for _, a in factors) else (-1) ** len(factors)


def totient_and_mobius(limit: int) -> tuple[list[int], list[int]]:
    phi = list(range(limit + 1))
    mu = [1] * (limit + 1)
    composite = bytearray(limit + 1)
    primes: list[int] = []
    for n in range(2, limit + 1):
        if not composite[n]:
            primes.append(n)
            phi[n] = n - 1
            mu[n] = -1
        for p in primes:
            product = n * p
            if product > limit:
                break
            composite[product] = 1
            if n % p == 0:
                phi[product] = phi[n] * p
                mu[product] = 0
                break
            phi[product] = phi[n] * (p - 1)
            mu[product] = -mu[n]
    return phi, mu


def add_column(
    basis: tuple[tuple[int, tuple[int, ...]], ...],
    column: tuple[int, ...],
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    """Extend a row-echelon basis of column vectors over F_101."""
    vector = list(column)
    for pivot, row in basis:
        multiplier = vector[pivot]
        if multiplier:
            for j in range(pivot, len(vector)):
                vector[j] = (vector[j] - multiplier * row[j]) % MODULUS
    pivot = next((j for j, value in enumerate(vector) if value), None)
    if pivot is None:
        return basis
    inverse = pow(vector[pivot], -1, MODULUS)
    normalized = tuple(value * inverse % MODULUS for value in vector)
    return tuple(sorted(basis + ((pivot, normalized),)))


def check_threshold(
    h: int, phi: list[int], mu: list[int], *, stop_at_four: bool = True
) -> tuple[int, tuple[int, ...], int, int, int]:
    budget = 4 * h
    rows = tuple(range(1, h, 2))
    candidates = sorted(
        (
            2 if e == 1 else phi[e],
            e,
            tuple(mu[e // d] % MODULUS if e % d == 0 else 0 for d in rows),
        )
        for e in range(1, budget * budget + 1, 2)
        if phi[e] <= budget
    )

    best_nullity = 0
    best_support: tuple[int, ...] = ()
    visited = 0

    def visit(
        start: int,
        remaining: int,
        basis: tuple[tuple[int, tuple[int, ...]], ...],
        support: tuple[int, ...],
    ) -> bool:
        nonlocal best_nullity, best_support, visited
        visited += 1
        nullity = len(support) - len(basis)
        if nullity > best_nullity:
            best_nullity, best_support = nullity, support
        if stop_at_four and nullity >= 4:
            return True
        for index in range(start, len(candidates)):
            weight, e, column = candidates[index]
            if weight > remaining:
                break
            if visit(
                index + 1,
                remaining - weight,
                add_column(basis, column),
                support + (e,),
            ):
                return True
        return False

    visit(0, budget, (), ())
    return best_nullity, best_support, visited, len(candidates), max(
        e for _, e, _ in candidates
    )


def check_counterexample(orders: tuple[int, ...], expected_cost: int) -> None:
    """Check explicit integer kernel vectors and a full-rank low minor.

    Divisor closure makes every omitted matrix row identically zero.
    The low-by-low minor is triangular with diagonal one, so the exact
    rational rank is its size; no inference from numerical rank is used.
    """
    h = min(orders)
    support = tuple(sorted(set().union(*(set(divisors(n)) for n in orders))))
    support_set = set(support)
    assert all(e % 2 for e in support)
    assert all(set(divisors(e)) <= support_set for e in support)
    low = tuple(e for e in support if e < h)
    high = tuple(e for e in support if e >= h)
    assert high == tuple(sorted(orders))
    cost = 2 + sum(exact_phi(e) for e in support if e > 1)
    assert cost == expected_cost <= 4 * h
    low_matrix = [
        [exact_mu(e // d) if e % d == 0 else 0 for e in low] for d in low
    ]
    assert all(low_matrix[j][j] == 1 for j in range(len(low)))
    assert all(low_matrix[j][k] == 0 for j in range(len(low)) for k in range(j))
    for n in orders:
        # Rows outside support are zero by the independently checked closure.
        assert all(
            sum(exact_mu(e // d) for e in divisors(n) if e % d == 0) == 0
            for d in low
        )
        assert [int(n % e == 0) for e in high] == [int(n == e) for e in high]
    assert len(support) - len(low) == 4


def check_counterexamples() -> None:
    for k in range(15, 79, 2):
        orders = (
            (k - 2) * (k + 2) * (k + 4),
            k * k * (k + 4),
            k * (k + 2) ** 2,
            k * (k + 2) * (k + 4),
        )
        cost = 4 * k**3 + 15 * k**2 - 4 * k - 23
        assert cost - 4 * min(orders) == 77 - (k - 6) ** 2 < 0
        check_counterexample(orders, cost)
    primes = (181, 193, 199)
    assert all(factorization(p) == ((p, 1),) for p in primes + (179, 191, 197))
    product = 181 * 193 * 199
    orders = tuple((p - 2) * (product // p) for p in primes) + (product,)
    assert all(all(a == 1 for _, a in factorization(n)) for n in orders)
    check_counterexample(orders, 4 * product - 3 * sum(product // p for p in primes) + 1)
    assert orders == (6874853, 6879629, 6881801, 6951667)
    print(
        "PASS: 32 odd-k fixtures and one squarefree fixture have exact "
        "rational nullity 4 within D(E)<=4h; universal nullity<=3 is false.",
        flush=True,
    )


def check_five_row_allocations() -> None:
    """Check actual widths, base parity and pair contacts, not just nullity."""
    for k in range(19, 51, 2):
        orders = (
            (k - 2) * (k + 2) * (k + 4),
            k * k * (k + 4),
            k * (k + 2) ** 2,
            k * (k + 2) * (k + 4),
        )
        support = tuple(sorted(set().union(*(set(divisors(n)) for n in orders))))
        raw = [tuple(int(n % e == 0) for e in support) for n in orders]
        raw.append(tuple(sum(raw[j][col] for j in range(3)) for col in range(len(support))))
        minimum = tuple(min(row[col] for row in raw) for col in range(len(support)))
        rows = [tuple(x - y for x, y in zip(row, minimum)) for row in raw]
        widths = tuple(max(row[col] for row in rows) for col in range(len(support)))
        assert len(set(rows)) == 5
        assert [row[0] for row in rows] == [0, 0, 0, 0, 2]
        assert widths[0] == 2
        degree = sum(exact_phi(e) * width for e, width in zip(support, widths))
        assert degree == 4 * k**3 + 15 * k**2 - k - 20
        assert 4 * min(orders) - degree == k * k - 15 * k - 44 > 0
        pair_contacts = []
        for left, right in combinations(rows, 2):
            difference = tuple(x - y for x, y in zip(left, right))
            moments = [
                sum(
                    exact_mu(e // d) * coefficient
                    for e, coefficient in zip(support, difference)
                    if e % d == 0
                )
                for d in support
            ]
            contact = next(d for d, value in zip(support, moments) if value)
            assert 4 * contact > degree
            pair_contacts.append(contact)
        assert min(pair_contacts) == min(orders)
    print(
        "PASS: 16 five-row allocation fixtures have even base parity, "
        "the exact width degree, and strict endpoint contact for all 10 pairs.",
        flush=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--max-h",
        type=int,
        default=19,
        help="largest odd contact threshold to enumerate (default: 19)",
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="enumerate the audited odd thresholds through h=31",
    )
    parser.add_argument(
        "--counterexamples-only",
        action="store_true",
        help="check the exact nullity-four counterexamples without enumeration",
    )
    args = parser.parse_args()
    check_counterexamples()
    check_five_row_allocations()
    if args.counterexamples_only:
        return
    max_h = 31 if args.full else args.max_h
    if max_h < 3 or max_h % 2 == 0:
        parser.error("--max-h must be an odd integer at least 3")
    max_index = (4 * max_h) ** 2
    phi, mu = totient_and_mobius(max_index)
    assert phi[1] == mu[1] == 1

    for h in range(3, max_h + 1, 2):
        maximum, witness, count, candidate_count, largest_index = check_threshold(
            h, phi, mu
        )
        print(
            f"h={h:2d} candidates={candidate_count:2d} "
            f"largest_e={largest_index:3d} supports={count:8d} "
            f"max_nullity_F{MODULUS}={maximum} witness={witness}",
            flush=True,
        )
        if maximum >= 4:
            raise SystemExit(
                f"Modular nullity >=4 at h={h}; rational-rank audit required: "
                f"{witness}"
            )
    print(
        f"PASS: modular nullity <=3 for every checked odd h<= {max_h}; "
        "therefore rational nullity <=3 on these finite budgets."
    )


if __name__ == "__main__":
    main()
