"""Exact prime-power bookkeeping for the composite-block moat transfer.

The probability estimate is proved in the companion note. This checks its
new arithmetic input by two independent calculations of the least rational
integer multiplier, including blocks with repeated and distinct primes.
"""

from itertools import product
from math import gcd

from check_norm_descent_rational_phase_gap import (
    add, conj, gaussian_gcd, mul, neg, norm, power,
)


def gaussian_product(values):
    result = (1, 0)
    for value in values:
        result = mul(result, value)
    return result


def check_block_families():
    families = [
        [power((2, 1), 3), power((3, 2), 2), power((4, 1), 2)],
        [mul((2, 1), (3, 2)), mul((4, 1), (5, 2)),
         mul((6, 1), (5, 4)), mul((7, 2), (8, 3))],
    ]
    cases = 0
    directions = 0
    for blocks in families:
        norms = [norm(block) for block in blocks]
        for j, block in enumerate(blocks):
            assert norm(gaussian_gcd(block, conj(block))) == 1
            assert all(gcd(norms[j], other) == 1 for other in norms[:j])
        p = 1
        for n in norms:
            p *= n
        signed_products = [gaussian_product(
            [conj(block) if sign else block for block, sign in zip(blocks, signs)])
            for signs in product((0, 1), repeat=len(blocks))]
        for x in range(-9, 10):
            for y in range(-9, 10):
                if gcd(x, y) != 1:
                    continue
                v = (x, y)
                all_contents = 1
                for b in signed_products:
                    assert norm(b) == p
                    content = norm(gaussian_gcd(v, b))
                    predicted = p // content
                    # b divides a*v iff both coordinates of a*v*bar(b)
                    # are divisible by N(b). This gives the true minimum
                    # without using Gaussian gcd or prime factorization.
                    numerator = mul(v, conj(b))
                    actual = p // gcd(p, gcd(*numerator))
                    assert actual == predicted
                    assert all((actual * c) % p == 0 for c in numerator)
                    all_contents *= content
                    cases += 1
                # Exact exponential form of E log G_sigma(v) <= log|v|.
                assert all_contents <= norm(v) ** (len(signed_products) // 2)
                directions += 1
    return cases, directions


def check_residual_update_failure():
    for k in range(1, 17):
        a = 2210 * k
        h, q1, q2 = (a, 1), (2 * a + 1, -2), (4 * a + 1, -4)
        d, n1, n2 = map(norm, (h, q1, q2))
        assert gcd(d, n1) == gcd(d, n2) == gcd(n1, n2) == 1
        p1, p2 = mul(h, q1), mul(h, q2)
        assert p1 == (2 * (a * a + 1) + a, 1)
        assert p2 == (4 * (a * a + 1) + a, 1)
        assert norm(gaussian_gcd(p1, conj(p1))) == 1
        assert norm(gaussian_gcd(p2, conj(p2))) == 1
        z0 = conj(gaussian_product((h, q1, q2)))
        z1 = gaussian_product((h, q1, conj(q2)))
        z2 = gaussian_product((h, q2, conj(q1)))
        assert norm(z0) == norm(z1) == norm(z2) == d * n1 * n2
        delta1, delta2 = add(z1, neg(z0)), add(z2, neg(z0))
        assert delta1 == mul((0, 2), conj(q2))
        assert delta2 == mul((0, 2), conj(q1))
        difference = add(delta1, neg(delta2))
        assert difference == mul((0, 4), h)
        # Stronger than distinct residues modulo the whole composite block:
        # the updates differ modulo every Gaussian prime dividing q1.
        assert norm(gaussian_gcd(difference, q1)) == 1
        assert norm(gaussian_gcd(z0, z1)) == n2
        assert norm(gaussian_gcd(z0, z2)) == n1
        assert norm(delta1) // n2 == norm(delta2) // n1 == 4
    return 16


if __name__ == "__main__":
    cases, directions = check_block_families()
    update_cases = check_residual_update_failure()
    print(f"PASS: {cases} exact least-multiplier identities;")
    print(f"      {directions} averaged-content inequalities, two composite families.")
    print(f"PASS: {update_cases} equal-residual, different-update circle fixtures.")
