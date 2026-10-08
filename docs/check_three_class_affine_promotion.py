#!/usr/bin/env python3
"""Exact finite checks for pairwise phase separation and slope promotion.

This checker supports the three-original-class-per-v_2-bucket argument in
three_class_affine_uniform_count.md. It uses exact arithmetic in Q(phi),
phi^2=phi+1. It is not a proof of the general theorem or of the packet
projection lemma; those are handled in the accompanying proof notes.
"""
from __future__ import annotations

from itertools import combinations, product
from math import gcd

# Elements a+b*phi of Z[phi], phi=(1+sqrt(5))/2.
Q = tuple[int, int]
ZERO: Q = (0, 0)
ONE: Q = (1, 0)
PHI: Q = (0, 1)


def add(x: Q, y: Q) -> Q:
    return x[0] + y[0], x[1] + y[1]


def neg(x: Q) -> Q:
    return -x[0], -x[1]


def sub(x: Q, y: Q) -> Q:
    return add(x, neg(y))


def mul(x: Q, y: Q) -> Q:
    a, b = x
    c, d = y
    # phi^2 = phi + 1
    return a * c + b * d, a * d + b * c + b * d


def power(x: Q, n: int) -> Q:
    assert n >= 0
    out = ONE
    while n:
        if n & 1:
            out = mul(out, x)
        x = mul(x, x)
        n >>= 1
    return out


def phi_power(n: int) -> Q:
    """Return phi**n in Q(phi), including negative integer n."""
    if n >= 0:
        return power(PHI, n)
    # phi^-1 = phi - 1
    return power((-1, 1), -n)


def is_zero(x: Q) -> bool:
    return x == ZERO


def conjugate(x: Q) -> Q:
    """Quadratic conjugation phi -> 1-phi."""
    a, b = x
    return a + b, -b


def determinant(x: Q, y: Q) -> int:
    """Determinant of x,y as columns in the rational basis (1,phi)."""
    return x[0] * y[1] - x[1] * y[0]


def v_p(n: int, p: int) -> int:
    assert n > 0 and p > 1
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def odd_primes(n: int) -> list[int]:
    assert n > 0
    while n % 2 == 0:
        n //= 2
    out = []
    d = 3
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 2
    if n > 1:
        out.append(n)
    return out


def lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


def lcm_all(values: tuple[int, ...] | list[int]) -> int:
    out = 1
    for value in values:
        out = lcm(out, value)
    return out


def ratio_distinct(slopes: tuple[int, ...], intercepts: tuple[int, ...]) -> bool:
    """No pair (A,B) is proportional, i.e. B/A ratios are distinct."""
    return all(intercepts[i] * slopes[j] != intercepts[j] * slopes[i]
               for i, j in combinations(range(len(slopes)), 2))


def phase_sign(A: int, B: int, n_parity: int, k: int) -> int:
    """The sign epsilon in gamma(A,B)^k = epsilon*i*phi^(-Bk).

    A is even, B and k are odd. Since i^(Bk)(-1)^(A*n*k/2)
    = i^((An+B)k), this is +i or -i.
    """
    assert A % 2 == 0 and B % 2 and k % 2 and n_parity in (0, 1)
    residue = ((A * n_parity + B) * k) % 4
    assert residue in (1, 3)
    return 1 if residue == 1 else -1


def cleared_phase(A: int, B: int, n_parity: int, k: int) -> Q:
    """Field scalar after stripping common i from gamma(A,B)^k."""
    sign = phase_sign(A, B, n_parity, k)
    return (sign * phi_power(-B * k)[0],
            sign * phi_power(-B * k)[1])


def cleared_frequency_term(A: int, B: int, n_parity: int,
                           k: int, moment: int) -> Q:
    """Return A*S*gamma^k after stripping common i.

    The raw log coefficient at physical frequency t=A*k is
    -2*(S/k)*gamma^k. Multiplying by -t/2 gives A*S*gamma^k.
    """
    phase = cleared_phase(A, B, n_parity, k)
    return phase[0] * A * moment, phase[1] * A * moment


def active_at_frequency(slopes: tuple[int, ...], t: int) -> list[int]:
    return [A for A in slopes if t % A == 0 and (t // A) % 2 == 1]


def promotion_prime(A: int, others: tuple[int, ...]) -> int:
    """Choose p with v_p(A)<v_p(B) for some other current slope B."""
    candidates = []
    for B in others:
        for p in odd_primes(B):
            if v_p(A, p) < v_p(B, p):
                candidates.append(p)
    assert candidates, (A, others)
    return min(candidates)


def promote_to_common(slopes: tuple[int, ...]) -> tuple[int, int, int]:
    """Iterate the exact minimum-slope promotion rule.

    Returns (number_steps, original_lcm, final_common_slope). Distinct slope
    values represent merged groups; duplicates are merged after every step.
    """
    current = sorted(set(slopes))
    H = lcm_all(current)
    steps = 0
    while len(current) > 1:
        A = current[0]
        others = tuple(current[1:])
        p = promotion_prime(A, others)
        assert p * A <= H and H % (p * A) == 0
        nxt = sorted(set(current[1:] + [p * A]))
        assert lcm_all(nxt) == H
        assert sum(nxt) > sum(current) or len(nxt) < len(current)
        current = nxt
        steps += 1
        assert steps <= 3 * H  # finiteness guard for these small fixtures
    return steps, H, current[0]


def check_phi_arithmetic() -> int:
    checks = 0
    assert mul(PHI, PHI) == add(PHI, ONE)
    for a in range(-80, 81):
        pa = phi_power(a)
        for b in range(-80, 81):
            pb = phi_power(b)
            assert mul(pa, pb) == phi_power(a + b)
            assert mul(pa, phi_power(-a)) == ONE
            checks += 2
    for n in range(-160, 161):
        p = phi_power(n)
        assert mul(p, conjugate(p)) == ((-1) ** n, 0)
        checks += 1
    return checks


def check_pairwise_phase_separation() -> int:
    """Check exact cleared-coefficient independence with varied intercepts."""
    cases = 0
    slopes_families = (
        (6, 10, 14),
        (10, 14, 22),
        (14, 22, 26),
        (18, 30, 42),
        (2, 6, 10),
    )
    intercept_choices = (
        (-9, -5, -1), (-9, 1, 7), (-5, 3, 11),
        (-1, 5, 11), (1, 3, 7), (5, -7, 9),
        (7, -11, 3),
    )
    for slopes in slopes_families:
        for Bs in intercept_choices:
            if not ratio_distinct(slopes, Bs):
                continue
            # Every active pair is checked at all odd common frequencies up
            # through a fixed range. This is finite evidence; the proof uses
            # the nonzero exponent difference and irrationality of phi^m.
            max_t = 1200
            for t in range(1, max_t + 1):
                active = [(A, t // A) for A in slopes
                          if t % A == 0 and (t // A) % 2 == 1]
                assert len(active) <= 3
                for (A, k), (C, ell) in combinations(active, 2):
                    # B/A != B'/C ensures exponents -Bk and -B'ell differ.
                    for parity in (0, 1):
                        phase1 = cleared_phase(A, Bs[slopes.index(A)], parity, k)
                        phase2 = cleared_phase(C, Bs[slopes.index(C)], parity, ell)
                        # Nonzero determinant means these two exact Q(phi)
                        # columns are linearly independent over Q. This
                        # certifies arbitrary rational moments, not just the
                        # small integer samples below.
                        assert determinant(
                            (A * phase1[0], A * phase1[1]),
                            (C * phase2[0], C * phase2[1]),
                        ) != 0
                        for s1, s2 in product(range(-1, 2), repeat=2):
                            value = add(
                                cleared_frequency_term(A, Bs[slopes.index(A)],
                                                      parity, k, s1),
                                cleared_frequency_term(C, Bs[slopes.index(C)],
                                                      parity, ell, s2),
                            )
                            if is_zero(value):
                                assert s1 == 0 and s2 == 0, (
                                    slopes, Bs, t, parity, (A, k), (C, ell),
                                    s1, s2, value,
                                )
                            cases += 1
    return cases


def check_clearing_denominator() -> int:
    checks = 0
    # At t=30, slopes 6 and 10 contribute at k=5 and k=3. The raw log
    # coefficient is -2*S*gamma^k/k; multiplying by -t/2 yields A*S.
    for A, B, k, S in ((6, 1, 5, 2), (10, 3, 3, -1), (14, 5, 5, 3)):
        t = A * k
        # Avoid integer-division assumptions: compare as exact rationals.
        from fractions import Fraction
        cleared_multiplier_q = Fraction(-t, 2) * Fraction(-2 * S, k)
        assert cleared_multiplier_q == A * S
        # It is generally not S; this catches dropping the slope multiplier.
        if A != 1:
            assert cleared_multiplier_q != S
        assert t == A * k
        checks += 1
    return checks


def check_frequency_exclusion_and_promotion() -> int:
    """Exhaust small same-v2 triples, promotion paths, and LCM invariant."""
    checked_triples = 0
    odd_parts = tuple(range(1, 32, 2))
    common_factors = (2, 8, 32)
    for common in common_factors:
        for ms in combinations(odd_parts, 3):
            slopes = tuple(common * m for m in ms)
            H = lcm_all(slopes)
            colors = {j: A for j, A in enumerate(slopes)}
            # For each promoted minimum and chosen prime, p-free frequencies
            # leave at most two original slope colors. The selected minimum
            # is always present at t=A*k.
            cur = tuple(sorted(set(slopes)))
            while len(cur) > 1:
                A = cur[0]
                p = promotion_prime(A, cur[1:])
                for k in range(1, min(4 * H // A, 1001), 2):
                    if k % p:
                        t = A * k
                        active_colors = [color for color, slope in colors.items()
                                         if t % slope == 0
                                         and (t // slope) % 2 == 1]
                        assert any(colors[color] == A for color in active_colors)
                        assert len(active_colors) <= 2, (
                            colors, A, p, k, active_colors)
                for color, slope in colors.items():
                    if slope == A:
                        colors[color] = p * A
                nxt = sorted(set(cur[1:] + (p * A,)))
                assert lcm_all(nxt) == H
                cur = tuple(nxt)
                assert set(colors.values()) == set(cur)
            steps, original_h, final = promote_to_common(slopes)
            assert original_h == H == final
            assert steps >= 1
            checked_triples += 1

    # Include a varied-incomparable three-slope path with explicit expected
    # LCM preservation and final common slope 210.
    slopes = (6, 10, 14)
    steps, H, final = promote_to_common(slopes)
    assert (H, final) == (210, 210)
    assert steps >= 3
    return checked_triples


def main() -> None:
    arithmetic = check_phi_arithmetic()
    phase = check_pairwise_phase_separation()
    denom = check_clearing_denominator()
    triples = check_frequency_exclusion_and_promotion()
    print(f"Q(phi) arithmetic checks: {arithmetic}")
    print(f"exact two-color phase checks: {phase}")
    print(f"cleared-denominator examples: {denom}")
    print(f"same-v2 slope triples promoted to their original LCM: {triples}")
    print("finite checks passed; general irrationality/promotion arguments require proof")


if __name__ == "__main__":
    main()
