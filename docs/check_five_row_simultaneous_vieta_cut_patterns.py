"""Exact combinatorics for simultaneous five-row Vieta selected cores.

This checker only verifies cut incidence and content exponents. It does
not construct endpoint points or bound the correction integers E_i.
"""

ROWS = tuple(range(5))
CUTS = tuple(mask for mask in range(1, 1 << 5))


def size(mask):
    return bin(mask).count("1")


def selected(mask, i):
    s = size(mask)
    inside = bool(mask & (1 << i))
    return (not inside and s >= 2) or (inside and s >= 4)


def h(mask, i):
    s = size(mask)
    inside = bool(mask & (1 << i))
    return max(0, 2 * (s - 1) - 5) if inside else max(0, 2 * s - 4)


def check_incidence():
    for i in ROWS:
        assert sum(selected(mask, i) for mask in CUTS) == 16
        for j in ROWS:
            if i != j:
                assert sum(selected(mask, i) and bool(mask & (1 << j))
                           for mask in CUTS) == 11

    for mask in CUTS:
        s = size(mask)
        multiplicity = sum(selected(mask, i) for i in ROWS)
        expected = {1: 0, 2: 3, 3: 2, 4: 5, 5: 5}[s]
        assert multiplicity == expected
        if s == 1:
            assert not any(selected(mask, i) for i in ROWS)
        elif s in (4, 5):
            assert all(selected(mask, i) for i in ROWS)


def check_content_exponents():
    for i in ROWS:
        assert sum(h(mask, i) for mask in CUTS) == 19

    total_by_size = {s: 0 for s in range(1, 6)}
    for mask in CUTS:
        total_by_size[size(mask)] = sum(h(mask, i) for i in ROWS)
    assert total_by_size == {1: 0, 2: 0, 3: 4, 4: 8, 5: 15}


def check_explicit_patterns():
    for i in ROWS:
        residual = [mask for mask in CUTS
                    if selected(mask, i) and size(mask) in (2, 3)]
        assert len(residual) == 10
        assert sum(size(mask) == 2 for mask in residual) == 6
        assert sum(size(mask) == 3 for mask in residual) == 4

    selected_multiplicity = {s: [] for s in range(1, 6)}
    for mask in CUTS:
        selected_multiplicity[size(mask)].append(
            sum(selected(mask, i) for i in ROWS))
    assert {s: sorted(values) for s, values in selected_multiplicity.items()} == {
        1: [0] * 5,
        2: [3] * 10,
        3: [2] * 10,
        4: [5] * 5,
        5: [5],
    }


if __name__ == "__main__":
    check_incidence()
    check_content_exponents()
    check_explicit_patterns()
    print("PASS: simultaneous selected cores have 0/3/2/5/5 incidences by cut size.")
    print("PASS: each row has 16 selected cuts, each retained row shares 11.")
    print("PASS: content exponents sum to 19 per row and 0,0,4,8,15 jointly.")
    print("PASS: formal common size-4/5 core has six blocks; residual rows have 6+4 cuts.")
