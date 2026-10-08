"""Exact finite checks for extreme-cut reflection profile bounds."""

from itertools import product


def cuts(k):
    """Unoriented nontrivial cuts, canonicalized by complement."""
    full = (1 << k) - 1
    return [m for m in range(1, full) if m < (full ^ m)]


def threshold_mask(e, t):
    return sum((x >= t) << i for i, x in enumerate(e))


def cut_class(mask, k):
    full = (1 << k) - 1
    comp = full ^ mask
    return min(mask, comp)


def check_extreme_bound(k, max_w):
    for e in product(range(max_w + 1), repeat=k):
        if min(e) != 0 or max(e) != max_w:
            continue
        for a in range(k):
            for b in range(a + 1, k):
                q = abs(e[a] + e[b] - max_w)
                extreme = 0
                full = (1 << k) - 1
                target = (1 << a) | (1 << b)
                for t in range(1, max_w + 1):
                    m = threshold_mask(e, t)
                    if m == target or m == (full ^ target):
                        extreme += 1
                assert q >= extreme, (k, e, a, b, q, extreme)


def check_uniform_counts(k):
    for ell in range(3, k + 1):
        # Enumerate every global unoriented cut and restrict it to J={0,...,
        # ell-1}; this checks the multiplicities rather than only the closed
        # formulas.  The retained pair is (0,1).
        full_j = (1 << ell) - 1
        extreme_j = cut_class((1 << 0) | (1 << 1), ell)
        nonconstant = 0
        extreme = 0
        for m in cuts(k):
            local = sum(((m >> i) & 1) << i for i in range(ell))
            if local == 0 or local == full_j:
                continue
            nonconstant += 1
            if cut_class(local, ell) == extreme_j:
                extreme += 1
        assert nonconstant == 2 ** (k - ell) * (2 ** (ell - 1) - 1)
        assert extreme == 2 ** (k - ell)
        assert extreme * (2 ** (ell - 1) - 1) == nonconstant


def check_cancellation_fixture():
    e = (1, 1, 0, 2)
    masks = [threshold_mask(e, t) for t in (1, 2)]
    assert [cut_class(m, 4) for m in masks] == [4, 7]
    assert abs(e[0] + e[1] - 2) == 0
    assert sum((e[0] >= t) == (e[1] >= t) for t in (1, 2)) == 2


def check_sharp_formal_profile(k):
    """Build the fixed-pair sharpness allocation with unit layer weights."""
    assert k >= 3
    a, b = 0, 1
    rest = list(range(2, k))
    target = {m: 0 for m in cuts(k)}
    q = 0

    # Complementary nonempty subsets of the remaining rows.  A two-layer
    # (0,1,2)-allocation supplies both cut classes and has q_ab=0.
    for mask in range(1, 1 << len(rest)):
        comp = ((1 << len(rest)) - 1) ^ mask
        if not comp or mask >= comp:
            continue
        e = [0] * k
        for j, row in enumerate(rest):
            e[row] = 0 if mask & (1 << j) else 2
        e[a] = e[b] = 1
        assert abs(e[a] + e[b] - 2) == 0
        for t in (1, 2):
            target[cut_class(threshold_mask(e, t), k)] += 1

    # The sole agreeing class represented by all remaining rows.
    e = [0] * k
    e[a] = e[b] = 1
    target[cut_class(threshold_mask(e, 1), k)] += 1
    q += 1

    # Every separating class gets a one-layer allocation and q_ab=0.
    for m in cuts(k):
        if ((m >> a) & 1) != ((m >> b) & 1):
            e = [(m >> i) & 1 for i in range(k)]
            target[m] += 1
            assert abs(e[a] + e[b] - 1) == 0

    assert set(target.values()) == {1}
    assert q == 1


def main():
    for k in range(3, 8):
        check_extreme_bound(k, 4)
        check_uniform_counts(k)
        check_sharp_formal_profile(k)
    check_cancellation_fixture()
    print("extreme-cut profile checks passed")


if __name__ == "__main__":
    main()
