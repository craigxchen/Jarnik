#!/usr/bin/env python3
"""Exact bounded search for Boolean widths of integer simplex rows.

Rows have integer coordinates, coordinate sum one, and are drawn from
[-B,B]^n.  Every reported quantity is integer arithmetic: if W is the sum of
the widths over all Boolean x, then F=W/2**n.
"""

from __future__ import annotations

import argparse
import itertools
from fractions import Fraction
from typing import Iterable, Sequence


def rows_in_box(n: int, bound: int) -> list[tuple[int, ...]]:
    """All rows in [-bound,bound]^n whose coordinates sum to one."""
    out: list[tuple[int, ...]] = []
    for prefix in itertools.product(range(-bound, bound + 1), repeat=n - 1):
        last = 1 - sum(prefix)
        if -bound <= last <= bound:
            out.append(prefix + (last,))
    return out


def dot_values(row: Sequence[int], n: int) -> tuple[int, ...]:
    """Dot products row*x in binary-mask order."""
    vals = []
    for mask in range(1 << n):
        vals.append(sum(v for i, v in enumerate(row) if mask >> i & 1))
    return tuple(vals)


def width_sum(rows: Sequence[Sequence[int]], values: Sequence[Sequence[int]] | None = None) -> int:
    """Sum_x(max_i row_i*x-min_i row_i*x), exactly as an integer."""
    n = len(rows[0])
    vv = values if values is not None else [dot_values(r, n) for r in rows]
    return sum(max(v[j] for v in vv) - min(v[j] for v in vv) for j in range(1 << n))


def rank(rows: Sequence[Sequence[int]]) -> int:
    """Exact rational rank, using fraction-free integer elimination."""
    a = [list(map(int, r)) for r in rows]
    if not a:
        return 0
    rank_ = 0
    n = len(a[0])
    for col in range(n):
        pivot = next((i for i in range(rank_, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank_], a[pivot] = a[pivot], a[rank_]
        p = a[rank_][col]
        # Ordinary exact elimination is sufficient for the small bounded scans.
        for i in range(len(a)):
            if i == rank_ or not a[i][col]:
                continue
            q = a[i][col]
            for j in range(col, n):
                a[i][j] = a[i][j] * p - q * a[rank_][j]
        rank_ += 1
        if rank_ == len(a):
            break
    return rank_


def affine_rank(rows: Sequence[Sequence[int]]) -> int:
    """Exact dimension of the affine span of the supplied rows."""
    if not rows:
        return 0
    base = rows[0]
    return rank([[x - y for x, y in zip(r, base)] for r in rows[1:]])


def affine_dependent(rows: Sequence[Sequence[int]]) -> bool:
    return affine_rank(rows) < len(rows) - 1


def pair_distance(a: Sequence[int], b: Sequence[int]) -> int:
    """Positive mass of a-b (equal to half its l1 norm when sums agree)."""
    return sum(max(x - y, 0) for x, y in zip(a, b))


def binary_equal_weight_rows(n: int, q: int) -> list[tuple[int, ...]]:
    """All incidence vectors of q-subsets of [n]."""
    return [tuple(1 if i in S else 0 for i in range(n)) for S in itertools.combinations(range(n), q)]


def audit_distance_two_configurations(n: int = 4) -> dict:
    """Audit the local four-row patterns used in the short-width argument.

    Translation by a common column baseline preserves widths, so it is enough
    to inspect equal-weight binary rows.  The scan keeps four rows whose pair
    distances are at most two and at least one is two.  For each pair-distance
    multiset it records the number of configurations and the least W.
    """
    by_pattern: dict[tuple[int, ...], dict] = {}
    checked = 0
    for q in range(1, n):
        rows = binary_equal_weight_rows(n, q)
        for idxs in itertools.combinations(range(len(rows)), 4):
            chosen = tuple(rows[i] for i in idxs)
            distances = tuple(sorted(pair_distance(a, b) for a, b in itertools.combinations(chosen, 2)))
            if distances[-1] != 2:
                continue
            checked += 1
            w = width_sum(chosen)
            entry = by_pattern.setdefault(distances, {"count": 0, "min_W": None, "example": None})
            entry["count"] += 1
            if entry["min_W"] is None or w < entry["min_W"]:
                entry["min_W"] = w
                entry["example"] = chosen
    return {"n": n, "checked": checked, "patterns": by_pattern}


def fmt(rows: Sequence[Sequence[int]]) -> str:
    return "{" + ", ".join("(" + ",".join(map(str, r)) + ")" for r in rows) + "}"


def scan(n: int, bound: int, m_values: Iterable[int]) -> dict:
    rows = rows_in_box(n, bound)
    values = [dot_values(r, n) for r in rows]
    result = {"n": n, "bound": bound, "row_count": len(rows), "m": {}}
    for m in m_values:
        if m > len(rows):
            continue
        best_all: tuple[int, tuple[int, ...]] | None = None
        best_dep: tuple[int, tuple[int, ...]] | None = None
        checked = dependent = general_violations = dependent_violations = 0
        general_rhs = (1 << n) * ((1 << (m - 1)) - 1)
        general_scale = 1 << (m - 1)
        for idxs in itertools.combinations(range(len(rows)), m):
            checked += 1
            chosen = tuple(rows[i] for i in idxs)
            w = width_sum(chosen, [values[i] for i in idxs])
            if w * general_scale < general_rhs:
                general_violations += 1
            if best_all is None or w < best_all[0]:
                best_all = (w, idxs)
            if affine_dependent(chosen):
                dependent += 1
                if w < (1 << n):
                    dependent_violations += 1
                if best_dep is None or w < best_dep[0]:
                    best_dep = (w, idxs)
        entry = {
            "checked": checked,
            "dependent_count": dependent,
            "general_violations": general_violations,
            "dependent_violations": dependent_violations,
        }
        if best_all is not None:
            w, idxs = best_all
            entry["all"] = {"W": w, "F": str(Fraction(w, 1 << n)), "rows": [rows[i] for i in idxs]}
        if best_dep is not None:
            w, idxs = best_dep
            entry["dependent"] = {"W": w, "F": str(Fraction(w, 1 << n)), "rows": [rows[i] for i in idxs]}
        result["m"][m] = entry
    return result


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, nargs="+", default=[2, 3])
    p.add_argument("--bound", type=int, default=2)
    p.add_argument("--m", type=int, nargs="+", default=None)
    p.add_argument("--audit", action="store_true", help="also audit equal-weight binary distance-two patterns")
    args = p.parse_args()
    for n in args.n:
        ms = args.m if args.m is not None else range(2, n + 2)
        print(scan(n, args.bound, ms))
        if args.audit:
            print(audit_distance_two_configurations(n))


if __name__ == "__main__":
    main()
