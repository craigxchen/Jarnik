"""Exact small checks for the Gaussian least-radius lcm/gcd formula."""

from fractions import Fraction
from itertools import product
from math import log


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def mul(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def conj(a):
    return a[0], -a[1]


def norm(a):
    return a[0] * a[0] + a[1] * a[1]


def nearest(num, den):
    """Nearest integer to num/den, with an arbitrary tie choice."""
    if num >= 0:
        return (2 * num + den) // (2 * den)
    return -((-2 * num + den) // (2 * den))


def divmod_gaussian(a, b):
    den = norm(b)
    raw = mul(a, conj(b))
    q = nearest(raw[0], den), nearest(raw[1], den)
    return q, add(a, (-mul(q, b)[0], -mul(q, b)[1]))


def gcd_gaussian(a, b):
    while b != (0, 0):
        _, r = divmod_gaussian(a, b)
        a, b = b, r
    return a


def exact_div(a, b):
    q, r = divmod_gaussian(a, b)
    assert r == (0, 0), (a, b, q, r)
    return q


def lcm_gaussian(a, b):
    return mul(exact_div(a, gcd_gaussian(a, b)), b)


def lcm_all(values):
    out = (1, 0)
    for value in values:
        out = lcm_gaussian(out, value)
    return out


def gcd_all(values):
    out = values[0]
    for value in values[1:]:
        out = gcd_gaussian(out, value)
    return out


def qmul(a, b):
    return tuple(Fraction(x) for x in mul(a, b))


def qdiv(a, b):
    n = mul(a, conj(b))
    d = norm(b)
    return tuple(Fraction(x, d) for x in n)


def gaussian_prime_candidates():
    # Small split Gaussian primes, found by direct norm primality.
    def prime(n):
        return n >= 2 and all(n % p for p in range(2, int(n**0.5) + 1))

    out = []
    for a in range(1, 20):
        for b in range(1, a + 1):
            p = a * a + b * b
            if prime(p):
                out.append((a, b))
    return out


def factor_product(factors):
    out = (1, 0)
    for factor, exponent in factors:
        for _ in range(exponent):
            out = mul(out, factor)
    return out


def valuation_rows(rows, factors):
    """Known valuations at each selected split prime orientation."""
    out = []
    for row in rows:
        row_values = {}
        for key, factor in factors.items():
            exponent = 0
            value = row
            while value != (0, 0):
                q, r = divmod_gaussian(value, factor)
                if r != (0, 0):
                    break
                exponent += 1
                value = q
            row_values[key] = exponent
        out.append(row_values)
    return out


def check_profile_and_corrections():
    m = 3
    all_blocks = list(product((0, 1), repeat=m))
    nonempty = [b for b in all_blocks if any(b)]
    factors = gaussian_prime_candidates()
    assert len(factors) >= len(nonempty) + 3

    # One distinct oriented Gaussian prime per nonempty Boolean block.
    block_factor = {b: factors[n] for n, b in enumerate(nonempty)}
    rows0 = []
    row_factors0 = []
    for i in range(m):
        chosen = [(block_factor[b], 1) for b in nonempty if b[i]]
        rows0.append(factor_product(chosen))
        row_factors0.append(chosen)

    L0 = lcm_all(rows0)
    G0 = gcd_all(rows0)
    # Every nonempty non-full block appears in the lcm/gcd quotient exactly
    # once; the full block is the common gcd and cancels.
    expected = 1
    for b in nonempty:
        if b != (1,) * m:
            expected *= norm(block_factor[b])
    assert norm(L0) == norm(G0) * expected

    # Add arbitrary conjugate-primitive row corrections on fresh primes.
    correction_factors = factors[len(nonempty):len(nonempty) + 3]
    correction_exponents = (2, 1, 3)
    corrections = [
        factor_product([(conj(correction_factors[i]), correction_exponents[i])])
        for i in range(m)
    ]
    rows = [mul(row, correction) for row, correction in zip(rows0, corrections)]
    L = lcm_all(rows)
    G = gcd_all(rows)

    selected = dict((str(b), p) for b, p in block_factor.items())
    selected.update((f"K{i}", conj(correction_factors[i]))
                    for i in range(m))
    old_vals = valuation_rows(rows0, selected)
    new_vals = valuation_rows(rows, selected)
    correction_vals = valuation_rows(corrections, selected)
    weighted_change = 0.0
    weighted_budget = 0.0
    for key, factor in selected.items():
        old_range = max(v[key] for v in old_vals) - min(v[key] for v in old_vals)
        new_range = max(v[key] for v in new_vals) - min(v[key] for v in new_vals)
        correction_range = max(v[key] for v in correction_vals) - min(v[key] for v in correction_vals)
        assert abs(new_range - old_range) <= correction_range
        weighted_change += abs(new_range - old_range) * log(norm(factor)) / 2
        weighted_budget += sum(v[key] for v in correction_vals) * log(norm(factor)) / 2
    assert weighted_change <= weighted_budget + 1e-12
    assert abs(log(norm(L) / norm(G)) / 2 -
               log(norm(L0) / norm(G0)) / 2) <= weighted_budget + 1e-12

    return rows, L, G


def check_multiplier_and_anchor(rows, L, G):
    # The unanchored generator beta=conj(L)/G is a Gaussian rational.  It
    # realizes every phase as an integral Gaussian coordinate.
    beta = qdiv(conj(L), G)
    for row in rows:
        value = qmul(qmul(beta, row), qdiv((1, 0), conj(row)))
        assert all(x.denominator == 1 for x in value), (row, value)

    # Adjoining the relative anchor 1 changes the lcm/gcd data to (L,1),
    # so the least radius is |L|.  The explicit multiplier conj(L) realizes
    # it, and no nonzero Gaussian integer has smaller modulus.
    anchored = [(1, 0)] + rows
    La = lcm_all(anchored)
    Ga = gcd_all(anchored)
    assert norm(La) == norm(L)
    assert Ga in ((1, 0), (-1, 0), (0, 1), (0, -1))
    beta_anchor = conj(L)
    for row in rows:
        value = qmul(qmul(beta_anchor, row), qdiv((1, 0), conj(row)))
        assert all(x.denominator == 1 for x in value), (row, value)


def main():
    rows, L, G = check_profile_and_corrections()
    check_multiplier_and_anchor(rows, L, G)
    print("Least-radius lcm/gcd formula passes exact Gaussian checks.")
    print("Primewise correction range changes are bounded by sum_i log|K_i|.")
    print("Adding anchor 1 forces gcd=1 and changes the least radius from |L|/|G| to |L|.")


if __name__ == "__main__":
    main()
