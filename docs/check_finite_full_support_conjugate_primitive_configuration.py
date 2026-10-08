#!/usr/bin/env python3
"""Exact checks for the finite conjugate-primitive three-row fixture."""

from math import gcd, isqrt


ROWS = (1_469_178, 31_686, 153_548)
Y = (1, 1, 1)

# (name, norm, Gaussian divisor (real, imag), incident zero-based rows)
BLOCKS = (
    ("123", 109, (10, 3), (0, 1, 2)),
    ("12", 157, (-11, 6), (0, 1)),
    ("23", 13, (3, -2), (1, 2)),
    ("13", 85, (-9, 2), (0, 2)),
    ("1", 1_483_897, (1109, 504), (0,)),
    ("2", 4_513, (-47, -48), (1,)),
    ("3", 195_749, (-385, -218), (2,)),
)

ROW_NORM_FACTORS = (
    (5, 17, 89, 109, 157, 16_673),
    (13, 109, 157, 4_513),
    (5, 13, 17, 61, 109, 3209),
)


def is_prime_trial(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def gaussian_norm(z: tuple[int, int]) -> int:
    return z[0] * z[0] + z[1] * z[1]


def gaussian_div_exact(
    z: tuple[int, int], h: tuple[int, int]
) -> tuple[int, int]:
    """Return z/h, asserting exact divisibility in Z[i]."""
    x, y = z
    a, b = h
    n = gaussian_norm(h)
    real_num = x * a + y * b
    imag_num = y * a - x * b
    assert real_num % n == 0 and imag_num % n == 0
    return real_num // n, imag_num // n


def main() -> None:
    assert all(x > 0 and x % 2 == 0 for x in ROWS)
    assert all(y == 1 for y in Y)
    row_norms = tuple(x * x + 1 for x in ROWS)

    for n, factors in zip(row_norms, ROW_NORM_FACTORS):
        product = 1
        for p in factors:
            assert is_prime_trial(p)
            assert p % 4 == 1
            product *= p
        assert product == n

    block_norms = tuple(block[1] for block in BLOCKS)
    for i, n in enumerate(block_norms):
        assert all(gcd(n, other) == 1 for other in block_norms[i + 1 :])
    for _, n, h, incident_rows in BLOCKS:
        assert gaussian_norm(h) == n
        assert gcd(h[0], h[1]) == 1
        for row in incident_rows:
            gaussian_div_exact((ROWS[row], Y[row]), h)

    for row in range(3):
        quotient = (ROWS[row], Y[row])
        for _, _, h, incident_rows in BLOCKS[:4]:
            if row in incident_rows:
                quotient = gaussian_div_exact(quotient, h)
        expected = BLOCKS[4 + row][2]
        assert quotient == expected
        assert gaussian_norm(quotient) == BLOCKS[4 + row][1]

    # Exact row norm gcds equal the prescribed common-plus-pair blocks.
    assert gcd(row_norms[0], row_norms[1]) == 109 * 157 == 17_113
    assert gcd(row_norms[1], row_norms[2]) == 109 * 13 == 1_417
    assert gcd(row_norms[0], row_norms[2]) == 109 * 85 == 9_265

    d12 = ROWS[0] - ROWS[1]
    d23 = ROWS[1] - ROWS[2]
    d13 = ROWS[0] - ROWS[2]
    assert (d12, d23, d13) == (1_437_492, -121_862, 1_315_630)
    assert (d12 // 17_113, d23 // 1_417, d13 // 9_265) == (84, -86, 142)
    assert all(d != 0 for d in (d12, d23, d13))

    # Any full four-row Boolean refinement must split every old block's
    # rational-prime support into nonempty exclude/include-4 pieces.
    for prime_block in (109, 157, 13):
        assert is_prime_trial(prime_block)
        assert not any(1 < d < prime_block and prime_block % d == 0 for d in range(2, prime_block))

    print("finite conjugate-primitive fixture: all exact checks passed")


if __name__ == "__main__":
    main()
