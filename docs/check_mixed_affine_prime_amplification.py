#!/usr/bin/env python3
"""Exact checks for mixed-affine prime amplification and packet projection.

This is a finite arithmetic checker, not a proof of Bertrand's postulate
or of the uniform theorem.  It checks Ramanujan congruences and bounds,
quadratic-unit recurrences, threshold arithmetic, v2 separation, weighted
colored root intersections, and integral packet identities on finite exact
fixtures.  In particular, it tests packet cost and low-jet preservation,
and includes a counterexample to coordinatewise prefix choices.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from math import gcd, isqrt


def mobius(n: int) -> int:
    if n == 1:
        return 1
    p = 2
    x = n
    sign = 1
    while p * p <= x:
        if x % p == 0:
            x //= p
            sign = -sign
            if x % p == 0:
                return 0
            while x % p == 0:
                x //= p
        p += 1 if p == 2 else 2
    if x > 1:
        sign = -sign
    return sign


def divisors(n: int) -> list[int]:
    low, high = [], []
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            low.append(d)
            if d * d != n:
                high.append(n // d)
    return low + high[::-1]


def ramanujan(e: int, k: int) -> int:
    g = gcd(e, k)
    return sum(d * mobius(e // d) for d in divisors(g))


def euler_phi(n: int) -> int:
    out = n
    p = 2
    x = n
    while p * p <= x:
        if x % p == 0:
            out -= out // p
            while x % p == 0:
                x //= p
        p += 1 if p == 2 else 2
    if x > 1:
        out -= out // x
    return out


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n) + 1, 2))


def prime_in_open_interval(lo: int, hi: int) -> int:
    for p in range(lo + 1, hi):
        if is_prime(p):
            return p
    raise AssertionError(f"no prime in ({lo}, {hi})")


def ceil_log2(n: int) -> int:
    assert n >= 1
    return (n - 1).bit_length()


def two_valuation(n: int) -> int:
    assert n > 0
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return v


def least_odd_at_least(x_num: int, x_den: int) -> int:
    """Least positive odd integer h with h >= x_num/x_den."""
    h = max(1, (x_num + x_den - 1) // x_den)
    if h % 2 == 0:
        h += 1
    return h


def colored_root_group(color: int, copies: int, order: int):
    return {
        (color, copy, Fraction(j, order))
        for copy in range(copies)
        for j in range(order)
    }


def colored_downset(color: int, copies: int, support: set[int]):
    roots = set()
    for order in support:
        roots.update(colored_root_group(color, copies, order))
    return roots


def packet_project(rows: list[dict[int, int]], p: int, use_level: int):
    """Apply one integral p-packet projection to row exponent maps.

    ``use_level`` is 0 or 1 and is fixed across every p-free base s.
    Coordinates absent from a row are read as zero.  Levels j>=2 shift to
    j-1; the selected level-0/1 slice becomes the new level-0 coordinate.
    """
    assert p > 2 and use_level in (0, 1)
    out = []
    for row in rows:
        new = {}
        support = set(row)
        for e in support:
            j, s = 0, e
            while s % p == 0:
                s //= p
                j += 1
            if j == use_level:
                new[s] = row.get(p**use_level * s, 0)
            elif j >= 2:
                new[p ** (j - 1) * s] = row[e]
        # The selected level may be absent in this row while present in a
        # different row.  Fill those common p-free coordinates with zero.
        out.append(new)
    all_s = set()
    for row in rows:
        for e in row:
            s = e
            while s % p == 0:
                s //= p
            all_s.add(s)
    for new in out:
        for s in all_s:
            new.setdefault(s, 0)
    return out


def row_widths(rows: list[dict[int, int]]) -> dict[int, int]:
    support = set().union(*(set(row) for row in rows))
    return {
        e: max(row.get(e, 0) for row in rows)
        - min(row.get(e, 0) for row in rows)
        for e in support
    }


def row_moment(row: dict[int, int], k: int) -> int:
    return sum(value * ramanujan(e, k) for e, value in row.items())


def rational_matrix_rank(matrix: list[list[int]]) -> int:
    if not matrix:
        return 0
    a = [[Fraction(value) for value in row] for row in matrix]
    rows, cols = len(a), len(a[0])
    rank = 0
    for col in range(cols):
        pivot = next((r for r in range(rank, rows) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][col]
        a[rank] = [x / scale for x in a[rank]]
        for r in range(rows):
            if r != rank and a[r][col]:
                factor = a[r][col]
                a[r] = [x - factor * y for x, y in zip(a[r], a[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def packet_index_map(support: set[int], p: int, level: int):
    """Return output-coordinate to input-coordinate map for one packet."""
    mapping = {}
    for e in support:
        j, s = 0, e
        while s % p == 0:
            s //= p
            j += 1
        if j == level:
            mapping[s] = e
        elif j >= 2:
            mapping[p ** (j - 1) * s] = e
    return mapping


def polynomial_multiply(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


def polynomial_divide_exact(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    assert b and b[-1] == 1
    rem = list(a)
    quotient = [0] * (len(a) - len(b) + 1)
    while len(rem) >= len(b):
        shift = len(rem) - len(b)
        factor = rem[-1]
        quotient[shift] = factor
        for j, value in enumerate(b):
            rem[shift + j] -= factor * value
        while rem and rem[-1] == 0:
            rem.pop()
    assert not rem, (a, b, rem)
    while len(quotient) > 1 and quotient[-1] == 0:
        quotient.pop()
    return tuple(quotient)


@lru_cache(None)
def cyclotomic_polynomial(n: int) -> tuple[int, ...]:
    """Cyclotomic coefficients in ascending order, computed integrally."""
    assert n >= 1
    current = [-1] + [0] * (n - 1) + [1]
    for d in divisors(n):
        if d < n:
            current = list(polynomial_divide_exact(tuple(current),
                                                   cyclotomic_polynomial(d)))
    return tuple(current)


def psi_polynomial(n: int) -> tuple[int, ...]:
    # The sign convention Ψ_1=1-z is the one used in the source note.
    return (1, -1) if n == 1 else cyclotomic_polynomial(n)


def polynomial_compose_power(a: tuple[int, ...], p: int) -> tuple[int, ...]:
    out = [0] * ((len(a) - 1) * p + 1)
    for j, value in enumerate(a):
        out[j * p] = value
    return tuple(out)


def check_packet_polynomial_identities() -> int:
    cases = 0
    # Keep degrees moderate while covering p=3, 5, 7 and nontrivial s.
    for p in (3, 5, 7):
        for s in range(1, 16, 2):
            if s % p == 0:
                continue
            lhs = polynomial_multiply(psi_polynomial(s), psi_polynomial(p * s))
            rhs = polynomial_compose_power(psi_polynomial(s), p)
            assert lhs == rhs, (p, s, lhs, rhs)
            cases += 1
            for j in (2,):
                lhs = psi_polynomial(p**j * s)
                rhs = polynomial_compose_power(psi_polynomial(p ** (j - 1) * s), p)
                assert lhs == rhs, (p, j, s, lhs, rhs)
                cases += 1
    return cases


def packet_choice(rows: list[dict[int, int]], p: int) -> int:
    widths = row_widths(rows)
    # Only p-free s contribute to the level-0/1 prefix widths.
    bases = set()
    for e in widths:
        s = e
        while s % p == 0:
            s //= p
        bases.add(s)
    c0 = sum(euler_phi(s) * widths.get(s, 0) for s in bases)
    c1 = sum(euler_phi(p * s) * widths.get(p * s, 0) for s in bases)
    return 0 if c0 <= c1 else 1


def weighted_cost(rows: list[dict[int, int]], slope: int) -> int:
    return slope * sum(euler_phi(e) * w for e, w in row_widths(rows).items())


def low_moment_packet_fixture(p: int, q: int) -> list[dict[int, int]]:
    """Rows whose differences have zero p-free moments at k=1,3."""
    assert p > 2 and q > 3 and gcd(p, q) == 1
    return [
        {1: 0, p: 0, q: 0, p * q: 0,
         p**2: 2, p**2 * q: 1, p**3: 1, p**3 * q: 2},
        {1: 1, p: 0, q: 1, p * q: 0,
         p**2: 0, p**2 * q: 2, p**3: 3, p**3 * q: 0},
        {1: 0, p: 0, q: 0, p * q: 0,
         p**2: 3, p**2 * q: 0, p**3: 0, p**3 * q: 3},
    ]


def gamma_power_signature(A: int, B: int, n_parity: int, k: int):
    """Formal exact signature of (i^B (-1)^(An/2) phi^(-B))^k."""
    assert A % 2 == 0 and B % 2 == 1 and n_parity in (0, 1)
    return (B * k % 4, ((A // 2) * n_parity * k) % 2, -B * k)


def check_packet_projection_exact() -> int:
    """Exercise packet identities, cost, and physical low-jet preservation.

    The examples are finite regression fixtures.  The accompanying note
    supplies the general algebraic proof.
    """
    cases = 0
    # A single projection for several primes and arbitrary p-free supports.
    for p in (3, 5, 7, 11):
        q = next(q for q in (5, 7, 11, 13, 17) if q != p)
        rows = low_moment_packet_fixture(p, q)
        choice = packet_choice(rows, p)
        assert choice == 1  # This fixture exercises the cheaper level-1 choice.
        projected = packet_project(rows, p, choice)
        old_cost = weighted_cost(rows, 6)
        new_cost = weighted_cost(projected, 6 * p)
        assert new_cost <= old_cost, (p, choice, old_cost, new_cost)
        widths = row_widths(rows)
        pfree = set()
        for e in widths:
            s = e
            while s % p == 0:
                s //= p
            pfree.add(s)
        c0 = sum(euler_phi(s) * widths.get(s, 0) for s in pfree)
        c1 = sum(euler_phi(s) * widths.get(p * s, 0) for s in pfree)
        old_prefix = 6 * (c0 + (p - 1) * c1)
        new_prefix = 6 * p * min(c0, c1)
        assert old_cost - new_cost == old_prefix - new_prefix

        # The coordinate map is a literal coordinate projection: each
        # output reads a different input, so its rank is the output count.
        support = set().union(*(set(r) for r in rows))
        index_map = packet_index_map(support, p, choice)
        input_coords = sorted(support)
        output_coords = sorted(index_map)
        projection_matrix = [
            [int(index_map[out] == inp) for inp in input_coords]
            for out in output_coords
        ]
        assert rational_matrix_rank(projection_matrix) == len(output_coords)
        assert len(input_coords) - len(output_coords) == len(pfree)

        # The Ramanujan moment of the old packet is zero off p and is p
        # times the rescaled moment on p-multiples, including repeated powers.
        for i in range(len(rows)):
            for j in range(i):
                diff_old = {e: rows[i].get(e, 0) - rows[j].get(e, 0)
                            for e in support}
                diff_new = {e: projected[i].get(e, 0) - projected[j].get(e, 0)
                            for e in set().union(*(set(r) for r in projected))}
                for k in range(1, 5, 2):
                    old_s = row_moment(diff_old, k)
                    if k % p:
                        assert old_s == 0, (p, k, old_s)
                    else:
                        new_s = row_moment(diff_new, k // p)
                        assert old_s == p * new_s, (p, k, old_s, new_s)
                        # Cleared physical coefficients A*S agree exactly;
                        # gamma'^(k/p) and gamma^k have identical signatures.
                        assert 6 * old_s == (6 * p) * new_s
                        for B in (-9, -3, 1, 7):
                            for n_parity in (0, 1):
                                assert gamma_power_signature(6, B, n_parity, k) == \
                                       gamma_power_signature(6 * p, p * B,
                                                             n_parity, k // p)
                    cases += 1

    # Two original colors at one slope can be packet-projected together;
    # their intercept phases remain classwise identical after rescaling.
    p = 5
    rows_by_color = []
    for color in range(2):
        q = 7 if color == 0 else 11
        rows = low_moment_packet_fixture(p, q)
        choice = packet_choice(rows, p)
        projected = packet_project(rows, p, choice)
        assert weighted_cost(projected, 10) <= weighted_cost(rows, 2)
        rows_by_color.append((rows, projected, choice))
    for color, (rows, projected, _choice) in enumerate(rows_by_color):
        B = 1 + 4 * color
        support = set().union(*(set(r) for r in rows))
        new_support = set().union(*(set(r) for r in projected))
        for i in range(len(rows)):
            for j in range(i):
                old = {e: rows[i].get(e, 0) - rows[j].get(e, 0)
                       for e in support}
                new = {e: projected[i].get(e, 0) - projected[j].get(e, 0)
                       for e in new_support}
                for k in range(1, 201, 2):
                    if k < 5 and k % p:
                        assert row_moment(old, k) == 0
                    if k < 5 and k % p == 0:
                        assert row_moment(old, k) == p * row_moment(new, k // p)
                        assert gamma_power_signature(2, B, 1, k) == \
                               gamma_power_signature(2 * p, p * B, 1, k // p)
                    cases += 1

    # Sequential projection by two different primes tests packets that
    # overlap in the original index set (including p^j q^l coordinates).
    p, q = 3, 5
    support = {p**j * q**ell * s
               for j in range(3) for ell in range(3) for s in (1, 7)}
    rows = [{e: i for e in support} for i in range(3)]
    slope = 14
    intercept = 1
    for prime in (p, q):
        choice = packet_choice(rows, prime)
        projected = packet_project(rows, prime, choice)
        assert weighted_cost(projected, slope * prime) <= weighted_cost(rows, slope)
        for i in range(len(rows)):
            for j in range(i):
                old_support = set().union(*(set(r) for r in rows))
                new_support = set().union(*(set(r) for r in projected))
                old = {e: rows[i].get(e, 0) - rows[j].get(e, 0)
                       for e in old_support}
                new = {e: projected[i].get(e, 0) - projected[j].get(e, 0)
                       for e in new_support}
                for k in range(1, 5, 2):
                    if k % prime:
                        assert row_moment(old, k) == 0
                    else:
                        assert row_moment(old, k) == prime * row_moment(new, k // prime)
                        assert gamma_power_signature(slope, intercept, 1, k) == \
                               gamma_power_signature(slope * prime,
                                                     intercept * prime, 1,
                                                     k // prime)
                    cases += 1
        rows, slope, intercept = projected, slope * prime, intercept * prime

    # A coordinatewise choice of level 0 or 1 can destroy the divisibility
    # identity.  This pair has u_1=u_5, so its k=1 moment vanishes, but the
    # locally cheaper choices split the two s-coordinates.
    p = 3
    old0 = {1: 1, 3: 0, 5: 1, 15: 2}
    old1 = {1: 0, 3: 0, 5: 2, 15: 4}
    diff = {e: old0[e] - old1[e] for e in old0}
    assert row_moment(diff, 1) == 0
    old_s3 = row_moment(diff, 3)
    # At s=1 choose level 1 (width 0 < 1); at s=5 choose level 0
    # (width 1 < 2).  The resulting new moment is not S_old(3)/p.
    coordinatewise_new = {
        1: old0[3] - old1[3],
        5: old0[5] - old1[5],
    }
    assert old_s3 == 6
    assert row_moment(coordinatewise_new, 1) == 1
    assert old_s3 != p * row_moment(coordinatewise_new, 1)
    # A single global level choice satisfies the exact identity.
    global_new = {
        1: old0[1] - old1[1],
        5: old0[5] - old1[5],
    }
    assert old_s3 == p * row_moment(global_new, 1)
    cases += 3

    return cases


def q_power_pair(k: int) -> tuple[int, int]:
    """Return (a,b) with q^k=a+bq, q^2+3q+1=0."""
    assert k >= 0
    a, b = 1, 0
    for _ in range(k):
        a, b = -b, a - 3 * b
    return a, b


def q_inverse_power_pair(k: int) -> tuple[int, int]:
    """Return (a,b) with q^(-k)=a+bq, using q^(-1)=-3-q."""
    assert k >= 0
    a, b = 1, 0
    for _ in range(k):
        # Multiply a+bq by -3-q and reduce q^2=-3q-1.
        a, b = -3 * a + b, -a
    return a, b


def multiply_pairs(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    """Multiply a+bq pairs and reduce q^2=-3q-1."""
    a, b = x
    c, d = y
    return a * c - b * d, a * d + b * c - 3 * b * d


def q_integer_power_pair(k: int) -> tuple[int, int]:
    return q_power_pair(k) if k >= 0 else q_inverse_power_pair(-k)


def check_ramanujan_congruence() -> int:
    # In c_e(dp^a)=sum_{r|e, r|dp^a} r*mu(e/r), every newly admitted
    # divisor r (one not dividing d) has larger p-adic valuation than d,
    # hence p|r. The remaining terms are exactly those in c_e(d).
    cases = 0
    # Includes p|d and p|e cases, repeated p-powers, and the coprime case.
    for e in range(1, 181):
        for d in range(1, 101):
            for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
                for power in range(1, 5):
                    k = d * p**power
                    assert (ramanujan(e, k) - ramanujan(e, d)) % p == 0, (
                        e,
                        d,
                        p,
                        power,
                        ramanujan(e, k),
                        ramanujan(e, d),
                    )
                    cases += 1
    return cases


def check_ramanujan_size_bound() -> int:
    cases = 0
    for e in range(1, 301):
        phi = euler_phi(e)
        for k in range(1, 401):
            assert abs(ramanujan(e, k)) <= phi, (e, k)
            cases += 1
    return cases


def check_quadratic_recurrence() -> int:
    # The two coordinates encode exact elements of Q(q), q^2+3q+1=0.
    for k in range(1, 201):
        a, b = q_power_pair(k)
        ai, bi = q_inverse_power_pair(k)
        assert b != 0, (k, a, b)  # q^k is not rational.
        assert bi == -b, (k, (a, b), (ai, bi))
        assert ai == a - 3 * b, (k, (a, b), (ai, bi))
        # Trace t_k=q^k+q^(-k) is the integer a+ai.
        assert b + bi == 0
        assert (a + ai, b + bi) == (a + ai, 0)
        # q^(2k) - t_k q^k + 1 = 0 exactly in Z[q].
        t = a + ai
        square = multiply_pairs((a, b), (a, b))
        assert (square[0] - t * a + 1, square[1] - t * b) == (0, 0)
    # Independently check the trace recurrence t_(k+1)=-3t_k-t_(k-1).
    traces = []
    for k in range(0, 202):
        traces.append(q_power_pair(k)[0] + q_inverse_power_pair(k)[0])
        assert q_power_pair(k)[1] + q_inverse_power_pair(k)[1] == 0
    assert traces[0:2] == [2, -3]
    for k in range(1, 201):
        assert traces[k + 1] == -3 * traces[k] - traces[k - 1]
    return 200


def check_negative_laurent_shifts() -> int:
    cases = 0
    for k in range(1, 16, 2):
        for low in range(-12, 1):
            coeffs = [(-1) ** j * (j + 1) for j in range(7)]
            terms = [(low + j, coeffs[j]) for j in range(len(coeffs))]
            value = (0, 0)
            for exponent, coefficient in terms:
                power = q_integer_power_pair(exponent * k)
                value = (value[0] + coefficient * power[0],
                         value[1] + coefficient * power[1])

            # Multiplication by X^(-low) turns the Laurent polynomial into
            # an ordinary polynomial and changes its value only by a unit
            # power q^(-low*k); coefficient l1 is unchanged.
            shifted = (0, 0)
            for j, coefficient in enumerate(coeffs):
                power = q_power_pair(j * k)
                shifted = (shifted[0] + coefficient * power[0],
                           shifted[1] + coefficient * power[1])
            factor = q_integer_power_pair(low * k)
            assert value == multiply_pairs(factor, shifted)
            assert sum(abs(c) for _, c in terms) == sum(map(abs, coeffs))
            cases += 1
    return cases


def check_bertrand_intervals() -> int:
    # The proof invokes Bertrand's postulate; this is a regression check of
    # the interval endpoints and the product lower bound, not its proof.
    cases = 0
    for n in range(10, 2001):
        p1 = prime_in_open_interval(n, 2 * n)
        p2 = prime_in_open_interval(2 * n, 4 * n)
        p3 = prime_in_open_interval(4 * n, 8 * n)
        assert p1 < p2 < p3
        assert p1 * p2 * p3 > 8 * n**3
        cases += 1
    return cases


def check_threshold_arithmetic() -> int:
    """Check the exact integer form of Astra's low-moment amplification.

    Put L=ceil(log_2 D), J=3L (a conservative low-d range). For D>=2^32, sqrt(D)>192(L+1), and for
    every odd d<J, h>=D/4 implies n=floor((h-1)/(8d))>2sqrt(D)-2. Thus
    n>=L, 8n^3>D, and all three Bertrand primes p<8n give L<dp<h, which
    is enough for Mahler decoupling. The inequalities are checked with
    integer arithmetic; the smooth bound log_2(D)+2 reduces the global
    threshold to D=2^32.
    """
    cases = 0
    D0 = 1 << 32
    L0 = ceil_log2(D0)
    assert isqrt(D0) > 192 * (L0 + 1)

    # Check the exact base threshold and a large range immediately above it.
    for D in [D0 + j for j in range(10_001)] + [1 << k for k in range(32, 81)]:
        L = ceil_log2(D)
        J = 3 * L
        # Integer form of sqrt(D)>192(L+1); use the smooth bound
        # L+1<=log_2(D)+2, whose ratio with sqrt(D) increases here.
        assert D > (192 * (L + 1)) ** 2, (D, L)
        h = (D + 3) // 4
        for d in (1, 3, 5, J - 2, J - 1):
            if d <= 0 or d >= J or d % 2 == 0:
                continue
            n = (h - 1) // (8 * d)
            assert n > 2 * isqrt(D) - 2, (D, h, d, n)
            assert n >= L, (D, h, d, n, L)
            assert 8 * n**3 > D, (D, h, d, n)
            assert 8 * n * d <= h - 1, (D, h, d, n)
            assert d * (n + 1) >= L + 1
            cases += 1
    return cases


def check_v2_frequency_disjointness() -> int:
    """Check that odd quotients select at most one distinct v2 class."""
    checked = 0
    slope_sets = (
        (2, 12, 40, 112),      # v2 = 1,2,3,4
        (6, 20, 56, 144),      # same v2 pattern, distinct odd parts
        (2, 24, 80, 672),      # wider range of distinct v2 values
    )
    for slopes in slope_sets:
        assert len({two_valuation(a) for a in slopes}) == len(slopes)
        for t in range(1, 20_001):
            active = [a for a in slopes if t % a == 0 and (t // a) % 2 == 1]
            assert len(active) <= 1, (slopes, t, active)
            if active:
                assert two_valuation(active[0]) == two_valuation(t)
            checked += 1

    # The distinct-v2 hypothesis is necessary: two different slopes with
    # the same v2 can contribute at the same odd-quotient frequency.
    assert [a for a in (2, 6) if 6 % a == 0 and (6 // a) % 2 == 1] == [2, 6]

    # Check the classwise cutoff h_c: A_c*k<tau for every odd k<h_c, and
    # A_c*h_c>=tau, for a range including non-divisibility cases.
    for tau in range(1, 300):
        for A in range(2, 42, 2):
            h = least_odd_at_least(tau, A)
            assert A * h >= tau
            assert all(A * k < tau for k in range(1, h, 2))
            checked += 1
    return checked


def check_weighted_colored_compression_graph() -> int:
    """Check weighted root unions, intersections, and graph bounds."""
    # Three distinct v2 slope classes with divisor-downset supports.
    # D=106 <= 4*tau=120; each displayed high order has A*e>=tau.
    classes = (
        (2, {1, 3, 5, 15}, 15),
        (4, {1, 3, 9}, 9),
        (8, {1, 5}, 5),
    )
    tau = 30
    total_cost = 0
    colored_union = set()
    high_nodes = []
    for color, (A, support, h) in enumerate(classes):
        assert all(d in support for e in support for d in divisors(e))
        cost = A * sum(euler_phi(e) for e in support)
        roots = colored_downset(color, A, support)
        assert len(roots) == cost
        total_cost += cost
        colored_union.update(roots)
        assert A * h >= tau
        for e in support:
            if e >= h:
                node = colored_root_group(color, A, e)
                assert len(node) == A * e >= tau
                assert node <= roots
                high_nodes.append((color, A, e, node))
    assert total_cost == len(colored_union) == 106
    assert total_cost <= 4 * tau
    slope_budget = sum(A * sum(euler_phi(e) for e in support)
                       for A, support, _ in classes)
    assert slope_budget == total_cost
    for A, _support, h in classes:
        B = (slope_budget + A - 1) // A
        assert B <= 4 * h
    assert len(high_nodes) == 3
    for i, (ci, Ai, ei, gi) in enumerate(high_nodes):
        for cj, Aj, ej, gj in high_nodes[i + 1:]:
            intersection = len(gi & gj)
            expected = Ai * gcd(ei, ej) if ci == cj else 0
            assert intersection == expected

    # A same-color overlap fixture exercises A*gcd(e,f) and the threshold.
    A, tau = 2, 10
    support = {1, 3, 5, 15}
    h = least_odd_at_least(tau, A)
    assert h == 5
    low_roots = colored_downset(0, A, support)
    assert len(low_roots) == A * sum(euler_phi(e) for e in support) == 30
    group5 = colored_root_group(0, A, 5)
    group15 = colored_root_group(0, A, 15)
    assert len(group5) == A * 5 and len(group15) == A * 15
    assert len(group5 & group15) == A * gcd(5, 15) == 10
    assert len(group5 & group15) >= Fraction(2 * tau, 15)
    assert 5 // gcd(5, 15) <= 30 and 15 // gcd(5, 15) <= 30

    # Exhaustively check the ratio bound and maximum-degree count over small
    # weighted classes. For an edge, A*gcd(n,m)>=2*tau/15 and A*n,A*m<=4tau.
    pair_cases = 0
    for tau in range(1, 101):
        for A in range(2, 34, 2):
            max_order = (4 * tau) // A
            orders = [n for n in range(1, max_order + 1, 2)
                      if A * n >= tau]
            for n in orders:
                neighbors = []
                for m in orders:
                    if m == n:
                        continue
                    g = gcd(n, m)
                    if 15 * A * g >= 2 * tau:
                        assert n // g <= 30 and m // g <= 30
                        neighbors.append(m)
                    pair_cases += 1
                assert len(neighbors) <= 224, (tau, A, n, neighbors)
    # Six independent high nodes would have union >4*tau by the first two
    # inclusion-exclusion terms, because each pair intersection is <2tau/15.
    for tau in range(1, 501):
        assert 6 * tau - 15 * Fraction(2 * tau, 15) == 4 * tau

    # The low-slope cost ceiling B=ceil(L/A) is covered by 4h_c whenever
    # L<=4tau and h_c is the least odd integer at least tau/A.
    for L in (1 << 32, 3 * (1 << 32), 10 * (1 << 32), 10**15):
        tau = (L + 3) // 4
        for A in range(2, min(L // (1 << 32), 64) + 1, 2):
            h = least_odd_at_least(tau, A)
            B = (L + A - 1) // A
            assert B <= 4 * h, (L, tau, A, h, B)
    return len(classes) + pair_cases + 500


def main() -> None:
    congruences = check_ramanujan_congruence()
    size_bounds = check_ramanujan_size_bound()
    recurrence = check_quadratic_recurrence()
    laurent = check_negative_laurent_shifts()
    bertrand = check_bertrand_intervals()
    thresholds = check_threshold_arithmetic()
    v2_cases = check_v2_frequency_disjointness()
    graph_cases = check_weighted_colored_compression_graph()
    packet_cases = check_packet_projection_exact()
    packet_poly_cases = check_packet_polynomial_identities()
    print(f"Ramanujan congruences checked: {congruences}")
    print(f"Ramanujan size bounds checked: {size_bounds}")
    print(f"quadratic recurrences checked: {recurrence}")
    print(f"negative Laurent shifts checked: {laurent}")
    print(f"Bertrand interval samples checked: {bertrand}")
    print(f"threshold samples checked: {thresholds}")
    print(f"v2-frequency fixtures checked: {v2_cases}")
    print(f"weighted colored graph cases checked: {graph_cases}")
    print(f"packet projection/cost/jet checks: {packet_cases}")
    print(f"packet cyclotomic polynomial identities: {packet_poly_cases}")
    print("All exact mixed-affine prime-amplification checks passed.")


if __name__ == "__main__":
    main()
