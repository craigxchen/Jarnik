"""Exact checks for signed-block witness height and the orientation-lcm gap.

The inequalities for arbitrary blocks and the infinite family are proved
in the notes. These finite checks retain prime powers, conjugation and
correction factors overlapping the original block supports.
"""

from itertools import combinations, product
from math import gcd

from check_norm_descent_rational_phase_gap import (
    add, conj, gaussian_gcd, mul, neg, norm, power,
)


def gproduct(values):
    result = (1, 0)
    for value in values:
        result = mul(result, value)
    return result


def least_integer_multiplier(v, b):
    # b | a*v iff N(b) divides both coordinates of a*v*bar(b).
    real, imag = mul(v, conj(b))
    return norm(b) // gcd(norm(b), gcd(real, imag))


def check_collinear_witnesses():
    families = [
        [power((2, 1), 2), (3, 2), (4, 1)],
        [mul((2, 1), (3, 2)), mul((4, 1), (5, 2)),
         mul((6, 1), (5, 4))],
    ]
    count = primitive_count = 0
    directions = [(x, y) for x in range(-6, 7) for y in range(-6, 7)
                  if gcd(x, y) == 1]
    signs = list(product((False, True), repeat=3))
    for blocks in families:
        for j, h in enumerate(blocks):
            assert norm(gaussian_gcd(h, conj(h))) == 1
            assert all(gcd(norm(h), norm(k)) == 1 for k in blocks[:j])
        signed = {s: gproduct([conj(h) if bit else h
                              for bit, h in zip(s, blocks)]) for s in signs}
        p = norm(signed[signs[0]])
        for sigma in signs:
            alpha = signed[sigma]
            # Deliberately overlap the correction with an incident block.
            correction = conj(blocks[0]) if sigma[0] else blocks[0]
            actual = mul(alpha, correction)
            assert norm(gaussian_gcd(actual, conj(actual))) == 1
            assert gcd(*actual) == 1
            for tau in signs:
                beta = signed[tau]
                q = 1
                for j in range(3):
                    if sigma[j] != tau[j]:
                        q *= norm(blocks[j])
                assert norm(gaussian_gcd(alpha, beta)) == p // q
                assert least_integer_multiplier(actual, beta) == q
                primitive_count += 1
                for v in directions:
                    ga = norm(gaussian_gcd(v, alpha))
                    gb = norm(gaussian_gcd(v, beta))
                    assert ga * gb <= norm(v) * (p // q)
                    ma = least_integer_multiplier(v, alpha)
                    mb = least_integer_multiplier(v, beta)
                    assert max(ma * ma, mb * mb) * norm(v) >= p * q
                    count += 1
    return count, primitive_count


def check_orientation_gap():
    cases = 0
    for r in range(2, 22, 2):
        for s in range(1, 13):
            d = r * r + 1
            a = 2 * r * d * s
            xp, xm = a + r, a - r
            hp, hm = (xp, 1), (xm, 1)
            qp, qm = norm(hp), norm(hm)
            fp = 4 * r * r * d * s * s + 4 * r * r * s + 1
            fm = 4 * r * r * d * s * s - 4 * r * r * s + 1
            assert (qp, qm) == (d * fp, d * fm)
            assert gcd(fp, fm) == 1 and gcd(qp, qm) == d
            assert norm(gaussian_gcd(hp, hm)) == 1
            assert norm(gaussian_gcd(hp, conj(hp))) == 1
            assert norm(gaussian_gcd(hm, conj(hm))) == 1
            ordinary = qp * qm // d
            assert ordinary % qp == ordinary % qm == 0
            assert qm * qm >= ordinary
            assert norm(mul(hp, hm)) == d * ordinary
            z0 = conj(mul(hp, hm))
            zp, zm = mul(hp, conj(hm)), mul(hm, conj(hp))
            assert len({z0, zp, zm}) == 3
            assert norm(z0) == norm(zp) == norm(zm) == d * ordinary
            assert norm(gaussian_gcd(gaussian_gcd(z0, zp), zm)) == 1
            # Actual third-edge denominator supplies the missing orientation.
            assert norm(gaussian_gcd(zp, zm)) == 1
            assert add(zp, neg(zm)) != (0, 0)
            cases += 1
    return cases


def exact_quotient(z, divisor):
    real, imag = mul(z, conj(divisor))
    denominator = norm(divisor)
    assert real % denominator == imag % denominator == 0
    return real // denominator, imag // denominator


def check_raw_flips():
    blocks = [power((2, 1), 2), (3, 2), (4, 1)]
    count = 0
    for signs in product((False, True), repeat=3):
        oriented = [conj(h) if bit else h for h, bit in zip(blocks, signs)]
        u = mul(gproduct(oriented), oriented[0])
        assert norm(gaussian_gcd(u, conj(u))) == 1
        for mask in range(1, 8):
            f = gproduct([h for j, h in enumerate(oriented) if mask & (1 << j)])
            q = norm(f)
            v = mul(exact_quotient(u, f), conj(f))
            a, b = f
            assert q * v[1] == (a*a-b*b)*u[1]-2*a*b*u[0]
            assert (abs(v[1])+abs(u[1]))**2*q*q >= 4*(q-1)*u[0]**2
            count += 1
    return count


def check_primitive_pell_half_flips():
    u, b = 2, 1
    for _ in range(32):
        # Omit the initial solution, whose second block is a unit.
        u, b = 9*u+20*b, 4*u+9*b
        assert u*u-5*b*b == -1 and u % 2 == 0 and b % 2 == 1
        f, a = (u+2*b, b), (b, 2*b-u)
        actual, raw = mul(f, a), mul(conj(f), a)
        corrected = mul((2, 1), raw)
        assert actual == (2*u*b, 1)
        assert raw == (4*b*b, 1-2*b*b)
        assert corrected == (10*b*b-1, 2)
        relation = add(mul((2, 0), actual), mul((2, -1), conj(raw)))
        relation = add(relation, neg(corrected))
        relation = add(relation, mul((-2, 0), conj(actual)))
        assert relation == (0, 0)
        q, n = norm(f), norm(a)
        assert q == 10*b*b+4*u*b-1 and n == 10*b*b-4*u*b-1
        assert q > 1 and n > 1 and gcd(q, n) == 1
        assert norm(actual) == q*n and norm(corrected) == 5*q*n
        for z in (f, a, actual, corrected):
            assert norm(gaussian_gcd(z, conj(z))) == 1
        # A fixed, exact comparability bound, independent of the index.
        assert n < q < 19*n
    return 32


def check_signed_circuit_incidence():
    signs = list(product((0, 1), repeat=3))
    count = boundary = 0
    for size in (2, 3, 4):
        for rows in combinations(signs, size):
            column_counts = [sum(row[j] for row in rows) for j in range(3)]
            singleton = any(n in (1, size-1) for n in column_counts)
            if size < 4:
                assert singleton
            elif not singleton:
                for j in range(3):
                    relative = tuple(row[j] ^ rows[0][j] for row in rows)
                    assert relative in ((0, 0, 0, 0), (0, 0, 1, 1),
                                        (0, 1, 0, 1), (0, 1, 1, 0))
                boundary += 1
            count += 1
    return count, boundary


if __name__ == "__main__":
    count, primitive_count = check_collinear_witnesses()
    gaps = check_orientation_gap()
    raw = check_raw_flips()
    pell = check_primitive_pell_half_flips()
    circuits, boundary = check_signed_circuit_incidence()
    print(f"PASS: {count} collinear height and content inequalities;")
    print(f"      {primitive_count} primitive-row multiplier identities with overlapping corrections.")
    print(f"PASS: {gaps} unbounded orientation-gap family fixtures.")
    print(f"PASS: {raw} exact raw-flip bounds; {pell} primitive Pell half-flip fixtures.")
    print(f"PASS: {circuits} sign-incidence cases, including {boundary} four-term boundaries.")
