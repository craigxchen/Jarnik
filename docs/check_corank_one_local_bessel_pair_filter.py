"""Exact sign/kernel/source checks for the local corank-one pair filter.

This finite checker does not verify the Bessel or Roth theorem inputs, and
its literal Gaussian rows do not have controlled angular span.
"""

from itertools import combinations, product
from math import isqrt


def chi(v, x):
    return -1 if bin(v & x).count("1") % 2 else 1


def bent_flipped(x):
    bit = lambda j: (x >> j) & 1
    parity = bit(0) * bit(1) + bit(2) * bit(3) + bit(4) * bit(5)
    value = -1 if parity % 2 else 1
    return -value if x == 0 else value


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    for k in range(3, isqrt(n) + 1, 2):
        if n % k == 0:
            return False
    return True


def split_prime_list(count):
    primes = []
    candidate = 100_000_001
    while len(primes) < count:
        if candidate % 4 == 1 and is_prime(candidate):
            primes.append(candidate)
        candidate += 4
    assert primes[-1] < 100_010_000
    return primes


def gaussian_prime(p):
    for a in range(1, isqrt(p) + 1):
        b2 = p - a * a
        if b2 <= 0:
            break
        b = isqrt(b2)
        if b > 0 and b * b == b2:
            return a, b
    raise AssertionError(f"No Gaussian representation for split prime {p}")


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def norm(z):
    return z[0] * z[0] + z[1] * z[1]


def power(z, exponent):
    result = (1, 0)
    while exponent:
        if exponent % 2:
            result = mul(result, z)
        z = mul(z, z)
        exponent //= 2
    return result


def unit_power(k):
    return [(1, 0), (0, 1), (-1, 0), (0, -1)][k % 4]


def full_fair_lattice_audit():
    checked = 0
    for m in range(2, 9):
        cuts = [mask for mask in range(1, (1 << m) - 1) if mask & 1]
        def length(c):
            return sum(abs(sum(c[i] for i in range(m) if mask >> i & 1))
                       for mask in cuts)
        assert length((1, -1) + (0,) * (m - 2)) == 2 ** (m - 2)
        if m <= 6:
            for c in product((-1, 0, 1), repeat=m):
                if sum(c) or not any(c):
                    continue
                assert length(c) >= 2 ** (m - 2)
                checked += 1
    return checked


def main():
    m = 64
    f = [bent_flipped(x) for x in range(m)]
    hat = [sum(f[x] * chi(v, x) for x in range(m)) for v in range(m)]
    assert hat[0] == 6
    assert set(hat[1:]) == {6, -10}
    assert sum(value * value for value in hat) == m * m

    weights = [abs(-hat[v] // 2) for v in range(1, m)] + [32]
    signs = [[(1 if -hat[v] > 0 else -1) * chi(v, x)
              for v in range(1, m)] + [f[x]] for x in range(m)]
    assert all({signs[x][k] for x in range(m)} == {-1, 1}
               for k in range(m))
    assert set(weights[:-1]) == {3, 5}
    assert all(sum(a * s for a, s in zip(weights, row)) == 3
               for row in signs)
    assert all(sum((signs[x][j] - signs[0][j]) * weights[j]
                   for j in range(m)) == 0 for x in range(1, m))

    # The 63 old Walsh columns already span the full row-difference space.
    assert all(sum(signs[x][j] for x in range(m)) == 0
               for j in range(m - 1))
    assert all(sum(signs[x][j] * signs[x][h] for x in range(m))
               == 64 * (j == h)
               for j in range(m - 1) for h in range(m - 1))

    # Old cuts are distinct and balanced, hence every two are incomparable.
    same = next((j, h) for j, h in combinations(range(m - 1), 2)
                if weights[j] == weights[h])
    j, h = same
    assert weights[j] == weights[h]
    assert 8 * 2 < m
    sj = {x for x in range(m) if signs[x][j] == 1}
    sh = {x for x in range(m) if signs[x][h] == 1}
    assert len(sj) == len(sh) == 32 and sj != sh
    assert not sj <= sh and not sh <= sj

    # The exact integer row certificate has only the chosen two columns.
    lam = [signs[x][j] - signs[x][h] for x in range(m)]
    assert sum(lam) == 0
    signature = [sum(lam[x] * signs[x][k] for x in range(m))
                 for k in range(m)]
    assert signature == [64 if k == j else -64 if k == h else 0
                         for k in range(m)]
    assert sum(map(abs, lam)) == 64  # phase cost = 64/(4*32)=1/2

    # Exact abstract obtuse Gram: 63 Walsh weights 1, extra weight 1/2.
    assert all(sum(signs[x][k] * signs[y][k] for k in range(m - 1)) == -1
               for x, y in combinations(range(m), 2))
    assert all(2 * sum(signs[x][k] * signs[y][k]
                       for k in range(m - 1))
               + signs[x][-1] * signs[y][-1] < 0
               for x, y in combinations(range(m), 2))

    # Distinct literal split-prime blocks, arbitrary row units, common d.
    old_primes = split_prime_list(m - 1)
    primes = old_primes + [5]
    gaussians = [gaussian_prime(p) for p in primes]
    assert all(norm(z) == p for z, p in zip(gaussians, primes))
    d = (7, 11)
    n_source = norm(d)
    for p in primes:
        n_source *= p
    rows = []
    for x in range(m):
        z = mul(d, unit_power(3 * x + 1))
        for s, prime in zip(signs[x], gaussians):
            z = mul(z, prime if s == 1 else (prime[0], -prime[1]))
        assert norm(z) == n_source
        rows.append(z)
    assert len(set(rows)) == m

    pairs = 0
    for x, y in combinations(range(m), 2):
        pair_norm = 1
        for k, p in enumerate(primes):
            if signs[x][k] != signs[y][k]:
                pair_norm *= p
        # Exact endpoint-derived pair norm for C=1/2: P^2>256 R^2.
        assert pair_norm * pair_norm > 256 * n_source
        pairs += 1

    # All other literal prime valuations cancel in the signed row product.
    # At j and h the pi-exponents are +32 and -32 respectively; conjugate
    # exponents are their negatives. The common factor cancels exactly.
    assert [sum(lam[x] * ((signs[x][k] + 1) // 2) for x in range(m))
            for k in range(m)] == [32 if k == j else -32 if k == h else 0
                                    for k in range(m)]
    numerator = (1, 0)
    denominator = (1, 0)
    for x, exponent in enumerate(lam):
        if exponent > 0:
            numerator = mul(numerator, power(rows[x], exponent))
        elif exponent < 0:
            denominator = mul(denominator, power(rows[x], -exponent))
    conjugate = lambda z: (z[0], -z[1])
    rhs_numerator = mul(power(gaussians[j], 32),
                        power(conjugate(gaussians[h]), 32))
    rhs_denominator = mul(power(conjugate(gaussians[j]), 32),
                          power(gaussians[h], 32))
    u = unit_power(sum(lam[x] * (3 * x + 1) for x in range(m)))
    assert mul(numerator, rhs_denominator) == mul(
        denominator, mul(u, rhs_numerator))
    fair_cases = full_fair_lattice_audit()

    print(f"PASS: M=64 full-support corank-one kernel weights {sorted(set(weights))}; "
          f"one equal-weight incomparable pair, exact Gaussian signed "
          f"row-product identity; "
          f"{pairs} literal Gaussian pair norms with C=1/2 condition.")
    print("Source angles are uncontrolled; Bessel and Roth are theorem inputs.")
    print(f"Full-fair lattice: {fair_cases} exact short coefficient checks; "
          "sharp pair-vector equality for M=2 through8.")


if __name__ == "__main__":
    main()
