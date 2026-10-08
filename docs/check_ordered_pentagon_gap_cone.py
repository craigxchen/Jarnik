#!/usr/bin/env python3
"""Exact cut-incidence check for ordered_pentagon_gap_cone_audit.md."""

from collections import Counter
from itertools import combinations


def audit(m: int) -> None:
    # The side not containing row zero uniquely represents an unoriented cut.
    cuts = [set(s) for r in range(1, m) for s in combinations(range(1, m), r)]
    pentagons = list(combinations(range(m), 5))
    columns = []
    for pentagon in pentagons:
        outer = {pentagon[0], pentagon[-1]}
        inner = set(pentagon[1:-1])
        retained_nonconstant = [
            cut
            for cut in cuts
            if 0 < len(set(pentagon) & cut) < 5
        ]
        assert len(retained_nonconstant) == 15 * 2 ** (m - 5)
        assert len(cuts) - len(retained_nonconstant) == 2 ** (m - 5) - 1
        column = [
            int(
                (outer <= cut and inner.isdisjoint(cut))
                or (inner <= cut and outer.isdisjoint(cut))
            )
            for cut in cuts
        ]
        assert sum(column) == 2 ** (m - 5)
        columns.append(column)

    multiplicities = [sum(column[j] for column in columns) for j in range(len(cuts))]
    assert len(cuts) == 2 ** (m - 1) - 1
    expected_blind = 10 * m - 36
    assert sum(value == 0 for value in multiplicities) == expected_blind
    assert all(multiplicities[j] == 0 for j, cut in enumerate(cuts) if min(len(cut), m-len(cut)) == 1)
    distribution = dict(sorted(Counter(multiplicities).items())) if m <= 7 else "checked"
    print(
        f"m={m}: pentagons={len(pentagons)}, cuts={len(cuts)}, "
        f"incidence={2 ** (m - 5)}, blind={expected_blind}, "
        f"multiplicities={distribution}"
    )


for row_count in range(5, 13):
    audit(row_count)
print("PASS")
