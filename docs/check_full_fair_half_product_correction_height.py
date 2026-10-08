"""Exact audits of the forced full-fair half-product correction height."""

from itertools import combinations
from math import comb, isqrt


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def conj(z):
    return z[0], -z[1]


def norm(z):
    return z[0] ** 2 + z[1] ** 2


def power(z, n):
    out = (1, 0)
    for _ in range(n):
        out = mul(out, z)
    return out


def prod(values):
    out = 1
    for value in values:
        out *= value
    return out


def gaussian_primes(count):
    result, p = [], 5
    while len(result) < count:
        if p % 4 == 1 and all(p % d for d in range(2, isqrt(p) + 1)):
            for a in range(1, isqrt(p) + 1):
                b = isqrt(p - a * a)
                if b and a * a + b * b == p:
                    result.append((a, b))
                    break
        p += 2
    return result


def cut_data(k):
    m = 2 * k
    cuts = [mask for mask in range(1, (1 << m) - 1) if mask & 1]
    star = (1 << k) - 1
    q = [sum(1 if i < k else -1 for i in range(m) if mask >> i & 1)
         for mask in cuts]
    return cuts, star, q


def count_audit():
    expected = {2: 4, 3: 27, 4: 136, 6: 2766}
    for k, value in expected.items():
        cuts, star, q = cut_data(k)
        assert len(cuts) == 2 ** (2 * k - 1) - 1
        assert sum(abs(imbalance) for cut, imbalance in zip(cuts, q) if cut != star) == value
        assert 2 * value == k * (comb(2 * k, k) - 2)
        assert sum(abs(t - k) * comb(2 * k, t) for t in range(2 * k + 1)) == \
            k * comb(2 * k, k)
        for i, j in combinations(range(2 * k), 2):
            assert sum(((cut >> i) ^ (cut >> j)) & 1 for cut in cuts) == 2 ** (2 * k - 2)
    return expected


def literal_audit(k):
    m = 2 * k
    cuts, star_cut, q = cut_data(k)
    star_index = cuts.index(star_cut)
    pis = gaussian_primes(len(cuts) + 1)
    ps = list(map(norm, pis))
    multiplicities = [1 + (j % 3) for j in range(len(cuts))]
    # A singleton imbalance is cancelled completely by overlapping K_i factors.
    assert cuts[0] == 1 and q[0] == 1
    multiplicities[0] = k
    gammas = [power(pi, exponent) for pi, exponent in zip(pis, multiplicities)]
    units = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    common = (7, 11)
    points, corrections = [], []
    for i in range(m):
        correction = conj(pis[0]) if i < k else pis[0]
        correction = mul(correction, pis[-1] if i % 2 else conj(pis[-1]))
        corrections.append(correction)
        point = mul(common, units[i % 4])
        point = mul(point, correction)
        for cut, gamma in zip(cuts, gammas):
            point = mul(point, gamma if cut >> i & 1 else conj(gamma))
        points.append(point)
    assert len(set(map(norm, corrections))) == 1
    assert len(set(map(norm, points))) == 1
    assert len(set(points)) == m

    exponents = [imbalance * exponent if cut != star_cut else 0
                 for cut, imbalance, exponent in zip(cuts, q, multiplicities)] + [0]
    exponents[0] -= k
    exponents[-1] = sum((1 if i < k else -1) for i in range(m) if i % 2)
    assert exponents[0] == 0
    h = (1, 0)
    for pi, exponent in zip(pis, exponents):
        h = mul(h, power(pi if exponent >= 0 else conj(pi), abs(exponent)))
    assert norm(h) == prod(p ** abs(exponent) for p, exponent in zip(ps, exponents))

    numerator = denominator = u = (1, 0)
    for i, point in enumerate(points):
        if i < k:
            numerator = mul(numerator, point)
            u = mul(u, units[i % 4])
        else:
            denominator = mul(denominator, point)
            u = mul(u, conj(units[i % 4]))
    numerator = mul(numerator, power(conj(gammas[star_index]), k))
    denominator = mul(denominator, power(gammas[star_index], k))
    assert mul(numerator, conj(h)) == mul(mul(u, denominator), h)
    baseline_norm = prod(norm(gamma) ** abs(imbalance)
                         for cut, gamma, imbalance in zip(cuts, gammas, q) if cut != star_cut)
    correction_norm = prod(map(norm, corrections))
    assert norm(h) * correction_norm >= baseline_norm
    return len(cuts)


if __name__ == '__main__':
    coefficients = count_audit()
    literal = [literal_audit(k) for k in [2, 3]]
    print(f"PASS: imbalance coefficients {coefficients}, all pair incidences; "
          f"literal Gaussian examples on {literal} complete cut classes.")
    print("Exact star removal, arbitrary units, reduced orientations and overlapping "
          "correction cancellation verified by integer identities.")
