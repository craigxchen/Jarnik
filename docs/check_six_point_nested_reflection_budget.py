#!/usr/bin/env python3
"""Exact finite checks for the sharp nested reflection profile."""

from __future__ import annotations

from collections import defaultdict


def cut_key(mask: int, m: int) -> int:
    full = (1 << m) - 1
    return min(mask, full ^ mask)


def all_cuts(m: int) -> set[int]:
    return {cut_key(mask, m) for mask in range(1, (1 << m) - 1)}


def threshold_cuts(e: tuple[int, ...]) -> dict[int, int]:
    """Return unordered threshold cuts with their level-gap weights."""
    m = len(e)
    levels = sorted(set(e))
    answer: dict[int, int] = defaultdict(int)
    for low, high in zip(levels, levels[1:]):
        mask = sum(1 << k for k, value in enumerate(e) if value >= high)
        answer[cut_key(mask, m)] += high - low
    return dict(answer)


def replacement(e: tuple[int, ...], a: int, b: int, c: int) -> tuple[int, ...]:
    result = list(e)
    result[c] = e[a] + e[b] - e[c]
    return tuple(result)


def width(e: tuple[int, ...]) -> int:
    return max(e) - min(e)


def build_profile(m: int) -> list[tuple[str, tuple[int, ...]]]:
    a, b, c = 0, 1, 2
    z_mask = sum(1 << k for k in range(3, m))
    profile: list[tuple[str, tuple[int, ...]]] = []

    # The T=empty boundary column: one binary cut {a,b}.
    profile.append(("empty", tuple(int(k in (a, b)) for k in range(m))))

    # Every nonempty T subset Z supplies the nested pair T and {a,b} union T.
    submask = z_mask
    while submask:
        e = tuple(2 if (submask >> k) & 1 else 1 if k in (a, b) else 0
                  for k in range(m))
        profile.append((f"T={submask}", e))
        submask = (submask - 1) & z_mask

    covered: set[int] = set()
    for _, e in profile:
        covered.update(threshold_cuts(e))

    # Orient every remaining cut away from c. These cuts contain exactly one
    # of a,b and are supplied by ordinary binary columns.
    full = (1 << m) - 1
    for key in sorted(all_cuts(m) - covered):
        mask = key if not ((key >> c) & 1) else full ^ key
        assert not ((mask >> c) & 1)
        assert ((mask >> a) & 1) + ((mask >> b) & 1) == 1
        e = tuple((mask >> k) & 1 for k in range(m))
        profile.append((f"binary={key}", e))
    return profile


def aggregate(profile: list[tuple[str, tuple[int, ...]]]) -> dict[int, int]:
    total: dict[int, int] = defaultdict(int)
    for _, e in profile:
        for cut, gap in threshold_cuts(e).items():
            total[cut] += gap
    return dict(total)


def check(m: int) -> None:
    a, b, c = 0, 1, 2
    z_mask = sum(1 << k for k in range(3, m))
    profile = build_profile(m)
    expected = {cut: 1 for cut in all_cuts(m)}
    assert aggregate(profile) == expected

    changed = [(name, width(replacement(e, a, b, c)) - width(e))
               for name, e in profile]
    assert sum(cost for _, cost in changed) == 0
    assert [(name, cost) for name, cost in changed if cost] == [
        ("empty", 1), (f"T={z_mask}", -1)
    ]

    new_profile = [(name, replacement(e, a, b, c)) for name, e in profile]
    assert aggregate(new_profile) == expected

    # Proper nested columns are permuted by T -> Z\T. The two boundary
    # columns trade the cut {a,b} against the pair {c},Z.
    by_name = dict(profile)
    for t_mask in range(1, 1 << m):
        if t_mask & ~z_mask or t_mask == z_mask:
            continue
        complement = z_mask ^ t_mask
        assert threshold_cuts(replacement(by_name[f"T={t_mask}"], a, b, c)) == \
               threshold_cuts(by_name[f"T={complement}"])
    assert set(threshold_cuts(replacement(by_name["empty"], a, b, c))) == {
        cut_key(1 << c, m), cut_key(z_mask, m)
    }
    assert set(threshold_cuts(replacement(by_name[f"T={z_mask}"], a, b, c))) == {
        cut_key((1 << a) | (1 << b), m)
    }

    print(f"PASS m={m}: {len(expected)} source and replacement cut weights are 1; "
          "total width change is 0.")


def main() -> None:
    for m in range(4, 9):
        check(m)


if __name__ == "__main__":
    main()
