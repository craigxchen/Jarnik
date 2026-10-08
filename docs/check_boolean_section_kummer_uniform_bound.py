#!/usr/bin/env python3
"""Exact arithmetic checks for the Boolean-section Kummer bound.

This is a finite regression checker for the closed formulas in the proof.
It is not the proof of the all-m bound; the monotone inequalities for m>=36
are elementary and are stated in the accompanying note.
"""

from math import comb


def section_counts(m, d):
    assert m >= 4 and d >= 1
    e = (2 ** (m - 1) - 2) * d
    t = {r: 2 * d * comb(m, r) for r in range(2, m)}
    f0 = sum(t.values())
    f1 = sum(r * t[r] for r in t)
    return e, t, f0, f1


def direct_k4(m, e, t):
    return m * (m - 4) * e - 8 * (m - 4) - sum(
        (r - 3) ** 2 * t[r] for r in t if r >= 3
    )


def closed_k4(m, d):
    return 32 - 8 * m + d * (
        (7 * m - 36) * 2 ** (m - 1) + m * m + 3 * m + 36
    )


def euler4(m, f1, t2):
    return 16 - 4 * m + f1 - t2


def closed_deficit(m, d):
    return 16 - 4 * m + d * (
        (36 - m) * 2 ** (m - 1) - 4 * m * m - 12 * m - 36
    )


def local_cluster(r, a):
    assert r >= 2 and a >= 1
    b = a % 2
    euler = 2 * a + (r - 2) * b
    cost = ((a + 1) // 2) * (r - 3) ** 2 + (a // 2) * (r - 1) ** 2
    kterm = (3 * r - 5) * a + (2 * r - 4) * b
    defterm = (11 - 3 * r) * a + (r - 2) * b
    return b, euler, cost, kterm, defterm


def repeated_closed_k_lower(m, d):
    return 32 - 8 * m + d * (
        (3 * m - 20) * 2 ** (m - 1) + 4 * m + 20
    )


def repeated_closed_def_upper(m, d):
    return 16 - 4 * m + d * (
        (36 - m) * 2 ** (m - 1) - 16 * m - 36
    )


def check_local_cluster_arithmetic():
    checked = 0
    for r in range(2, 41):
        for a in range(1, 21):
            b, euler, cost, kterm, defterm = local_cluster(r, a)
            assert kterm + cost == a * r * (r - 1)
            assert kterm + defterm == 3 * euler
            assert defterm <= (9 - 2 * r) * a
            assert kterm >= (3 * r - 5) * a
            checked += 1
    return checked


def sample_partitions(total):
    choices = [
        (total,),
        (1,) * total,
        (total - 1, 1),
    ]
    if total % 2 == 0:
        choices.append((total // 2, total // 2))
    return choices


def profile_totals(m, d, partitions_by_r):
    """Aggregate exact local data over all subset-size clusters."""
    assert m >= 4 and d >= 1
    total_a_pair = 0
    euler_local = 0
    cost_total = 0
    kterm_total = 0
    defterm_total = 0
    pair_count = 0
    for r in range(2, m):
        parts = partitions_by_r[r]
        assert sum(parts) == 2 * d
        ways = comb(m, r)
        a_pair = sum(a * comb(r, 2) for a in parts)
        pair_count += ways * a_pair
        total_a_pair += ways * sum(a * r * (r - 1) for a in parts)
        euler_local += ways * sum(local_cluster(r, a)[1] for a in parts)
        cost_total += ways * sum(local_cluster(r, a)[2] for a in parts)
        kterm_total += ways * sum(local_cluster(r, a)[3] for a in parts)
        defterm_total += ways * sum(local_cluster(r, a)[4] for a in parts)
    e = (2 ** (m - 1) - 2) * d
    assert pair_count == comb(m, 2) * e
    assert total_a_pair == m * (m - 1) * e

    # These are the exact global Euler, canonical-square, Kummer, and
    # deficit identities after the weighted pair-count cancellation.
    euler4_value = 16 - 4 * m + euler_local
    k4_value = 32 - 8 * m - 3 * m * e + kterm_total
    assert k4_value == m * (m - 4) * e - 8 * (m - 4) - cost_total
    def_value = 3 * euler4_value - k4_value
    expected_def = 16 - 4 * m + 3 * m * e + defterm_total
    assert def_value == expected_def

    # Equivalent K expression via local Kterms and the pair identity.
    assert k4_value == kterm_total - 3 * m * e - 8 * (m - 4)
    return e, euler4_value, k4_value, def_value, cost_total


def check_repeated_profiles():
    checked = 0
    for m in list(range(4, 13)) + [35, 36, 37, 64]:
        for d in range(1, 6):
            choices = sample_partitions(2 * d)
            # Exercise merged, simple, near-merged, and balanced partitions;
            # vary the choice with r to ensure aggregation handles mixed data.
            for shift in range(len(choices)):
                profile = {
                    r: choices[(r + shift) % len(choices)]
                    for r in range(2, m)
                }
                e, _, k4, deficit, cost = profile_totals(m, d, profile)
                max_cost = sum(
                    comb(m, r) * local_cluster(r, 2 * d)[2]
                    for r in range(2, m)
                )
                assert cost <= max_cost
                assert k4 >= m * (m - 4) * e - 8 * (m - 4) - max_cost
                assert k4 >= repeated_closed_k_lower(m, d)

                def_bound = 16 - 4 * m + 3 * m * e + sum(
                    comb(m, r) * 2 * d * (9 - 2 * r)
                    for r in range(2, m)
                )
                assert deficit <= def_bound
                assert def_bound == repeated_closed_def_upper(m, d)
                if m >= 36:
                    assert k4 >= repeated_closed_k_lower(m, d) > 0
                    assert deficit <= repeated_closed_def_upper(m, d) < 0
                checked += 1
    return checked


def check_base_genus_shifts():
    checked = 0
    genera = (0, 1, 2, 10)
    for m in (4, 7, 12, 35, 36, 37, 64):
        for d in (1, 2, 5):
            choices = sample_partitions(2 * d)
            profile = {r: choices[r % len(choices)] for r in range(2, m)}
            e, _, _, _, _ = profile_totals(m, d, profile)
            euler_local = sum(
                comb(m, r)
                * sum(local_cluster(r, a)[1] for a in profile[r])
                for r in range(2, m)
            )
            kterm_total = sum(
                comb(m, r)
                * sum(local_cluster(r, a)[3] for a in profile[r])
                for r in range(2, m)
            )
            defterm_total = sum(
                comb(m, r)
                * sum(local_cluster(r, a)[4] for a in profile[r])
                for r in range(2, m)
            )
            for g in genera:
                euler_g = 4 * (m - 4) * (g - 1) + euler_local
                k_g = 8 * (m - 4) * (g - 1) - 3 * m * e + kterm_total
                deficit_g = 3 * euler_g - k_g
                assert deficit_g == 4 * (m - 4) * (g - 1) + 3 * m * e + defterm_total
                if m >= 36:
                    k_lower = 8 * (m - 4) * (g - 1) + repeated_closed_k_lower(m, d)
                    def_upper = (
                        4 * (m - 4) * (g - 1)
                        - d * ((m - 36) * 2 ** (m - 1) + 16 * m + 36)
                    )
                    assert k_g >= k_lower
                    assert deficit_g <= def_upper
                checked += 1
    return checked


def check_formulas():
    checked = 0
    for m in range(4, 81):
        for d in range(1, 6):
            e, t, f0, f1 = section_counts(m, d)
            # Every section label contributes choose(r,2) pair incidences.
            assert sum(comb(r, 2) * t[r] for r in t) == comb(m, 2) * e
            assert f0 == sum(2 * d * comb(m, r) for r in range(2, m))
            assert f1 == sum(2 * d * r * comb(m, r) for r in range(2, m))

            k4a = direct_k4(m, e, t)
            k4b = closed_k4(m, d)
            assert k4a == k4b, (m, d, k4a, k4b)

            chi4 = euler4(m, f1, t[2])
            deficit = 3 * chi4 - k4a
            assert deficit == closed_deficit(m, d), (m, d, deficit)

            if m >= 36:
                assert k4a > 0, (m, d, k4a)
                assert deficit < 0, (m, d, deficit)
                assert repeated_closed_k_lower(m, d) > 0
                assert repeated_closed_def_upper(m, d) < 0
            checked += 1
    return checked


def check_subset_labels(max_m=10):
    """Check local inertia labels in F_2^(m-1) for the Kummer covers."""
    count = 0
    for m in range(4, max_m + 1):
        labels = [1 << j for j in range(m - 1)] + [(1 << (m - 1)) - 1]
        assert len(set(labels)) == m
        assert all(label != 0 for label in labels)
        xor_all = 0
        for label in labels:
            xor_all ^= label
        assert xor_all == 0

        # The sole nontrivial dependency is the sum of all m labels.
        subset_xors = {}
        for mask in range(1, 1 << m):
            value = 0
            members = []
            for j, label in enumerate(labels):
                if (mask >> j) & 1:
                    value ^= label
                    members.append(j)
            subset_xors[mask] = value
            if mask != (1 << m) - 1:
                assert value != 0, (m, mask)

        # A proper nonempty subset S of size at least two, together with
        # any member i in S, gives two nonzero distinct inertia labels.
        for mask, value in subset_xors.items():
            size = bin(mask).count("1")
            if 2 <= size <= m - 1:
                assert value != 0
                for j, label in enumerate(labels):
                    if (mask >> j) & 1:
                        assert value != label, (m, mask, j)

        # The original double-cover labels are themselves nonzero and
        # pairwise distinct, as needed for each pair C_i,C_j.
        for i in range(m):
            for j in range(i + 1, m):
                assert labels[i] != 0 and labels[j] != 0
                assert labels[i] != labels[j]
        count += 2**m - 2
    return count


def main():
    checked = check_formulas()
    labels = check_subset_labels()
    local = check_local_cluster_arithmetic()
    repeated = check_repeated_profiles()
    genus = check_base_genus_shifts()
    print(f"PASS: {checked} exact (m,d) formula cases; m=4..80, d=1..5.")
    print(f"PASS: F2^(m-1) local inertia labels and subset tests through m=10 ({labels} subsets).")
    print(f"PASS: {local} local repeated-root cluster cases; r=2..40, a=1..20.")
    print(f"PASS: {repeated} mixed subset-size multiplicity profiles; exact Euler/K/deficit and global bounds.")
    print(f"PASS: {genus} base-genus shifts for g=0,1,2,10 on selected profiles.")
    print("The finite checker verifies identities and bookkeeping, not the all-m proof.")


if __name__ == "__main__":
    main()
