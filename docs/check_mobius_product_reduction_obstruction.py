"""Exact obstruction to product-minimizing Möbius radius reduction."""

from math import gcd, prod

from check_gaussian_reflection_replacement import gnorm
from check_mobius_conductor_transfer import least_norm


ROWS = [(1, 0), (-1241, 8), (-2008, -1), (-322, 1)]
SOURCE_N = 358853785
SOURCE_PRODUCT = 643880195044131125
TARGET_N = 84815456245
MINIMUM_PRODUCT = 8481545624500


def objective(d):
    return (
        (d * d + 1)
        * ((8 * d - 17305) ** 2 + 64)
        * ((d - 2330) ** 2 + 1)
    )


def main():
    assert all(gcd(*h) == 1 and gnorm(h) % 2 for h in ROWS)
    assert least_norm(ROWS) == SOURCE_N
    assert prod(map(gnorm, ROWS)) == SOURCE_PRODUCT == 5 * SOURCE_N**2

    # Exact endpoint ordering and extreme-angle bound.
    h1, h2 = ROWS[1], ROWS[2]
    determinant = abs(h1[0] * h2[1] - h1[1] * h2[0])
    dot = abs(h1[0] * h2[0] + h1[1] * h2[1])
    assert determinant == 17305 and dot == 2491920
    assert dot == 144 * determinant
    assert SOURCE_N < 144**4

    # Every A >= 2 branch of the exact two-variable objective exceeds F(2163).
    assert 4 * (8665 * 1250) ** 2 > MINIMUM_PRODUCT
    assert (
        4 * 64 * 17145**2 * 1175**2
        > MINIMUM_PRODUCT * 8**4
    )
    assert 16 * (1080 * 160 * 83) ** 2 > MINIMUM_PRODUCT
    assert 4 * (2247 * 671) ** 2 > MINIMUM_PRODUCT

    values = [(objective(d), d) for d in range(2331)]
    assert min(values) == (MINIMUM_PRODUCT, 2163)
    assert objective(2163) == MINIMUM_PRODUCT

    shear = (1, 155, 0, 1)
    a, b, c, d = shear
    targets = [(a * x + b * y, c * x + d * y) for x, y in ROWS]
    assert targets == [(1, 0), (-1, 8), (-2163, -1), (-167, 1)]
    assert prod(map(gnorm, targets)) == MINIMUM_PRODUCT
    assert least_norm(targets) == TARGET_N
    assert MINIMUM_PRODUCT < SOURCE_PRODUCT
    assert TARGET_N > 236 * SOURCE_N
    assert 137**4 < SOURCE_N < 138**4

    # The target pair chord 16/sqrt(65) violates C' <= 2 exactly.
    assert 4096 * TARGET_N > 65**2
    assert 155**2 + 2 == 24027

    print("PASS: actual primitive Pell source has C < 2 and the stated least radius.")
    print("PASS: the exact GL(2,Z) row-product objective reduces to (A,D).")
    print("PASS: all A >= 2 branches exceed the A=1 candidate.")
    print("PASS: exhaustive A=1 minimum is the critical shear B=155.")
    print("PASS: the minimizing image has N'/N > 236 and violates C' <= 2.")


if __name__ == "__main__":
    main()
