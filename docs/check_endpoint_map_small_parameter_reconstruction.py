"""Exact source-only recovery of triangular endpoint-map candidates.

The algebraic recovery is tested without assuming endpoint hypotheses;
the small bounds themselves follow from the cited endpoint theorems.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt, lcm, prod

from check_gaussian_reflection_replacement import gconj, gmul, gnorm
from check_mobius_conductor_transfer import fixtures, least_norm
from check_mobius_reciprocal_stretch_grid import clips, gram, realize


def primitive(z):
    c = gcd(*z)
    return z[0] // c, z[1] // c


def reanchor(rows, j):
    return [primitive(gmul(h, gconj(rows[j]))) for h in rows]


def recover(anchor, E, Z, a):
    """Integer-only integrality and square tests; None rejects a candidate."""
    r, s = gmul(Z, gconj(anchor))
    if r % E or s % E:
        return None
    r, s = r // E, s // E
    if s % (2 * a):
        return None
    b = s // (2 * a)
    d2 = a * a - b * b - r
    if d2 <= 0:
        return None
    d = isqrt(d2)
    if d * d != d2 or gcd(a, b, d) != 1:
        return None
    return a, b, d


def all_edge_norm(rows):
    out = 1
    for (x, y), (u, v) in combinations(rows, 2):
        D, F = x * u + y * v, x * v - y * u
        c = gcd(D, F)
        D, F = D // c, F // c
        epsilon = 2 if D % 2 and F % 2 else 1
        out = lcm(out, (D * D + F * F) // epsilon)
    return out


def source_recovery_audit():
    sources = fixtures() + [[(1, 0), (1, 1), (2, 1), (3, 1)],
                           [(1, 0), (-26, 7), (3, -1), (-17, -1)]]
    matrices = [(a, b, d) for a in range(1, 4) for d in range(1, 4)
                for b in range(-4, 5)
                if gcd(a, b, d) == 1 and (a != d or b)]
    cases = clipped_cases = nontrivial_E = descents = increases = 0
    inverse_cases = 0
    for rows in sources:
        source, N = realize(rows)
        for j, anchor in enumerate(source):
            anchored_rows = reanchor(rows, j)
            assert anchored_rows[j] == (1, 0)
            assert least_norm(anchored_rows) == N
            for a, b, d in matrices:
                T, W = gram(a, b, d)
                Q = prod(p ** e for p, e in clips(a, b, d, anchor, N).items())
                C = gmul(W, anchor)
                assert N % Q == C[0] % Q == C[1] % Q == 0
                E, Z = N // Q, (C[0] // Q, C[1] // Q)
                assert Z != (0, 0)
                assert recover(anchor, E, Z, a) == (a, b, d)
                assert E * E * gnorm(W) == N * gnorm(Z)
                Fminus, Fplus = E * (T - 2 * a * d), E * (T + 2 * a * d)
                assert 0 < Fminus < Fplus
                assert Fminus * Fplus == N * gnorm(Z)
                assert Fplus - Fminus == 4 * E * a * d
                # The exact norm inequality used for the Z size bound.
                assert gnorm(Z) * N < T * T * E * E
                images = [(a * x + b * y, d * y) for x, y in anchored_rows]
                target, Nprime = realize(images)
                assert all_edge_norm(images) == Nprime == least_norm(images)
                inverse_Q = prod(p ** e for p, e in
                                 clips(d, -b, a, target[j], Nprime).items())
                assert inverse_Q == Q
                _, inverse_W = gram(d, -b, a)
                inverse_C = gmul(inverse_W, target[j])
                inverse_Z = inverse_C[0] // Q, inverse_C[1] // Q
                assert inverse_C[0] % Q == inverse_C[1] % Q == 0
                assert recover(target[j], Nprime // Q, inverse_Z, d) == (d, -b, a)
                cases += 1
                inverse_cases += 1
                clipped_cases += Q > 1
                nontrivial_E += E > 1
                descents += Nprime < N
                increases += Nprime > N
    assert descents and increases and clipped_cases and nontrivial_E
    return cases, inverse_cases, clipped_cases, descents, increases


def rejection_and_search_audit():
    """A finite exact search checks reconstruction in both directions."""
    anchor, N = (2, 1), 5
    accepted = 0
    recovered = set()
    for E in (1, 5):
        Q = N // E
        for x in range(-50, 51):
            for y in range(-50, 51):
                if (x, y) == (0, 0):
                    continue
                Z = x, y
                for a in range(1, 4):
                    candidate = recover(anchor, E, Z, a)
                    if candidate is None:
                        continue
                    _, b, d = candidate
                    _, W = gram(a, b, d)
                    assert W != (0, 0)
                    assert gmul(W, anchor) == (Q * x, Q * y)
                    recovered.add(candidate)
                    accepted += 1
    # E=N (Q=1) always supplies the unconditional algebraic representation.
    # Small-parameter endpoint bounds, not recovery alone, make the theorem useful.
    for a in range(1, 4):
        for b in range(-3, 4):
            for d in range(1, 4):
                if gcd(a, b, d) != 1 or (a == d and b == 0):
                    continue
                _, W = gram(a, b, d)
                Z = gmul(W, anchor)
                assert max(map(abs, Z)) <= 50
                assert (a, b, d) in recovered
    assert recover(anchor, 5, (1, 0), 1) is None  # nonintegral W
    assert recover(anchor, 1, (1, 0), 1) is None  # nonintegral b
    return accepted, len(recovered)


def constants(m):
    k, F, Qm = m // 4, (m - 1) ** 2 // 4, m * m // 4
    q = Fraction(m * (m - 1), 2 * k * ((m - 1) * (m - 2 * k - 1) + 1))
    g = Fraction(m, 4 * F)
    ell = Fraction(1, 4 * (m - 1)) + Fraction(Qm, m * (m - 1)) * q
    h = ell + 2 * g
    log2D = 2 + Fraction(m * (m - 1), F) + m * (ell + g)
    L = 2 * ((m + 1) // 2)
    u = Fraction(1, 2) + (Fraction(1, L) if L % 4 == 2
                          else Fraction(L, L * L - 4))
    rho = 1 - Fraction(4, m)
    return q, h, u, rho, log2D


def exponent_audit():
    for m in range(8, 10001):
        q, h, u, rho, log2D = constants(m)
        assert 0 < q < 1
        Aexponent = h / rho
        Zexponent = q - Fraction(1, 2) + (u + 2 * h) / rho
        v = q + Aexponent + 2 * Zexponent
        assert v == 3 * q - 1 + (2 * u + 5 * h) / rho
        if m >= 256:
            assert q <= Fraction(8, m)
            assert h <= Fraction(7, m)
            assert u <= Fraction(1, 2) + Fraction(2, m)
            assert log2D <= 14
            assert Aexponent <= Fraction(8, m)
            assert Zexponent <= Fraction(27, m)
            assert log2D + 4 * h / rho <= 15
            # Bound log2(5) by 3 in the constant for Zmax.
            assert 3 + u + m * q + 2 * log2D + 4 * (u + 2 * h) / rho <= 43
            assert v <= Fraction(70, m)
    q, h, u, rho, _ = constants(100000)
    v = 3 * q - 1 + (2 * u + 5 * h) / rho
    assert abs(100000 * v - Fraction(137, 4)) < Fraction(1, 100)
    multiplier_exponent = 2 * q - 1 + 2 * (u + 2 * h) / rho
    gap_exponent = q + 2 * h / rho
    assert abs(100000 * multiplier_exponent - 27) < Fraction(1, 100)
    assert abs(100000 * gap_exponent - Fraction(21, 2)) < Fraction(1, 100)
    return 10001 - 8


def main():
    cases, inverses, clips_count, descents, increases = source_recovery_audit()
    accepted, distinct = rejection_and_search_audit()
    exponents = exponent_audit()
    print(f"PASS: {cases} actual-source-anchor recoveries and {inverses} inverse recoveries;")
    print(f"      {clips_count} nontrivial clipped conductors, {descents} strict radius descents,")
    print(f"      {increases} strict increases, all with exact full image radii.")
    print(f"PASS: small search accepted {accepted} records for {distinct} distinct maps.")
    print(f"PASS: {exponents} rational constant/exponent checks and asymptotic coefficient 137/4.")


if __name__ == '__main__':
    main()
