"""Primewise checks for all-anchor core, lcm, and correction identities.

Each case is one split prime in a nonempty endpoint cut on four labels
0,1,2,3. Independent prime supports make these local checks composable.
The proof is in all_anchor_rotation_radius_cocycle.md; finite checks
guard its signs and exact correction exponents.
"""

from itertools import product


def check_one(e: int, r: int, mask: int, allocations: tuple[int, ...]) -> None:
    """Check one prime of original exponent e and trimmed exponent e-2r."""
    core = e - 2 * r
    assert core > 0
    assert len(allocations) == 4
    assert mask in range(1, 8)

    a0 = allocations[0]
    signed = [0] + [a - a0 for a in allocations[1:]]
    inside = [False] + [bool(mask & (1 << (i - 1))) for i in range(1, 4)]
    residual = [signed[i] - core * int(inside[i]) for i in range(4)]

    # E_0=lcm(bar K_i): the two Gaussian orientations are the positive
    # and negative residual exponents, respectively.
    ebar = max(0, *residual)
    epi = max(0, *(-value for value in residual))
    low = min(allocations)
    high = max(allocations)

    for anchor, allocation in enumerate(allocations):
        # L_a=z_a/gcd(z_0,z_1,z_2,z_3), in the pi/bar(pi) orientations.
        actual = (allocation - low, high - allocation)
        core_point = (core, 0) if inside[anchor] else (0, core)
        correction = (epi + residual[anchor], ebar - residual[anchor])
        predicted = (
            core_point[0] + correction[0],
            core_point[1] + correction[1],
        )
        assert min(*correction) >= 0, (e, r, mask, allocations, anchor)
        assert predicted == actual, (e, r, mask, allocations, anchor)
        assert sum(actual) == core + epi + ebar
        assert epi + ebar <= 2 * r  # Norm(E) divides D^2.

    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            cut = int(inside[j]) - int(inside[i])
            if cut:
                # The pair's core numerator divides its primitive
                # numerator with the correct Gaussian orientation.
                assert (signed[j] - signed[i]) * cut >= core, (
                    e,
                    r,
                    mask,
                    allocations,
                    i,
                    j,
                )


def main() -> None:
    count = 0
    for e in range(2, 9):
        for r in range((e - 1) // 2 + 1):
            low = range(r + 1)
            high = range(e - r, e + 1)
            for mask in range(1, 8):
                choices = [high if mask & (1 << i) else low for i in range(3)]
                for allocations in product(low, *choices):
                    check_one(e, r, mask, allocations)
                    count += 1
    assert count == 6573
    print(f"verified {count} primewise endpoint allocation fixtures")


if __name__ == "__main__":
    main()
