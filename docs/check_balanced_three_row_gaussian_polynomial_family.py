"""Exact polynomial and progression certificates, using only Python stdlib."""
from itertools import combinations
from math import gcd, prod

FACTORS = {'g': (0, -1), '12': (-6, -5), '13': (-6, 3),
           '23': (-6, 4), '1': (-4, 3), '2': (-3, 2), '3': (-7, -6)}
INCIDENT = [('g', '12', '13', '1'), ('g', '12', '23', '2'),
            ('g', '13', '23', '3')]
CONTENTS = [15, 2, 30]
CORRECTIONS = [(7, 16), (47, 90), (5, 14)]
S_CERT = {2, 3, 5, 7, 13, 17, 29, 37, 41, 53, 61, 101}
MODULUS = 8088851149545216141000


def factor(n):
    ans = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            ans[d] = ans.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        ans[n] = ans.get(n, 0) + 1
    return ans


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def norm(z):
    return z[0] * z[0] + z[1] * z[1]


def conjugate(z):
    return z[0], -z[1]


def gaussian_gcd(a, b):
    while b != (0, 0):
        denominator = norm(b)
        numerator = mul(a, conjugate(b))
        quotient = tuple((2 * x + denominator) // (2 * denominator) for x in numerator)
        subtract = mul(quotient, b)
        a, b = b, (a[0] - subtract[0], a[1] - subtract[1])
    return a


def exact_divide(a, b):
    x, y = mul(a, conjugate(b))
    denominator = norm(b)
    assert x % denominator == 0 and y % denominator == 0
    return x // denominator, y // denominator


def gaussian_lcm(a, b):
    return exact_divide(mul(a, b), gaussian_gcd(a, b))


def determinant(z, w):
    return z[0] * w[1] - w[0] * z[1]


def gaussian_product(items):
    out = (1, 0)
    for z in items:
        out = mul(out, z)
    return out


def polynomial_product(a, b):
    out = [(0, 0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            u, v = mul(x, y)
            q, r = out[i + j]
            out[i + j] = q + u, r + v
    return out


def polynomial_eval(coefficients, t):
    real = sum(z[0] * t ** i for i, z in enumerate(coefficients))
    imag = sum(z[1] * t ** i for i, z in enumerate(coefficients))
    return real, imag


def primitive(z):
    return gcd(abs(z[0]), abs(z[1])) == 1 and norm(z) % 2 == 1


def main():
    rows = []
    for blocks in INCIDENT:
        row = [(1, 0)]
        for label in blocks:
            row = polynomial_product(row, [FACTORS[label], (1, 0)])
        rows.append(row)
    expected_real = [[105, -256, 106, -16, 1],
                     [94, -195, 95, -15, 1],
                     [150, -439, 151, -19, 1]]
    assert [[z[0] for z in row] for row in rows] == expected_real
    assert [[z[1] for z in row] for row in rows] == [
        [240, 0, 0, 0, 0], [180, 0, 0, 0, 0], [420, 0, 0, 0, 0]]
    constants = {label: norm(z) for label, z in FACTORS.items()}
    quads = {label: [norm(z), 2 * z[0], 1] for label, z in FACTORS.items()}
    for (i, j), label, coefficient in [((0, 1), '12', -60),
                                       ((0, 2), '13', 180),
                                       ((1, 2), '23', 240)]:
        det_poly = [rows[i][d][0] * rows[j][0][1]
                    - rows[j][d][0] * rows[i][0][1] for d in range(5)]
        expected = polynomial_product([(x, 0) for x in quads['g']],
                                      [(x, 0) for x in quads[label]])
        assert det_poly == [coefficient * z[0] for z in expected]
    for i in range(3):
        assert gcd(*rows[i][0]) == CONTENTS[i]
        assert tuple(c // CONTENTS[i] for c in rows[i][0]) == CORRECTIONS[i]
        assert primitive(CORRECTIONS[i])
    assert [norm(z) for z in CORRECTIONS] == [305, 10309, 221]
    correction_gcd_norms = [norm(gaussian_gcd(CORRECTIONS[i], CORRECTIONS[j]))
                            for i, j in combinations(range(3), 2)]
    assert correction_gcd_norms == [61, 1, 13]
    assert norm(gaussian_gcd(gaussian_gcd(*CORRECTIONS[:2]), CORRECTIONS[2])) == 1
    ss = {2}
    for a, b in FACTORS.values():
        ss.update(factor(abs(b)))
        ss.update(factor(a * a + b * b))
    for z in CORRECTIONS:
        ss.update(factor(norm(z)))
    for (c, b, _), (e, d, _) in combinations(quads.values(), 2):
        u, v = d - b, e - c
        resultant = c * u * u - b * u * v + v * v
        assert resultant != 0
        ss.update(factor(abs(resultant)))
    assert ss == S_CERT
    modulus = prod(p ** (1 + max(factor(c).get(p, 0) for c in constants.values()))
                   for p in ss)
    assert modulus == MODULUS
    correction_lcm = (912, -211)
    circle_corrections = [correction_lcm, (-464, 813), (-348, 869), (-572, 741)]
    assert norm(correction_lcm) == 876265
    computed_lcm = (1, 0)
    for kk in CORRECTIONS:
        computed_lcm = gaussian_lcm(computed_lcm, conjugate(kk))
    assert norm(exact_divide(computed_lcm, correction_lcm)) == 1
    assert norm(correction_lcm) * prod(correction_gcd_norms) == prod(norm(z) for z in CORRECTIONS)
    assert prod(constants.values()) == 4500 * norm(correction_lcm)
    for i, kk in enumerate(CORRECTIONS):
        assert mul(circle_corrections[i + 1], conjugate(kk)) == mul(correction_lcm, kk)
        assert norm(circle_corrections[i + 1]) == norm(correction_lcm)
    for h in [1, 2, 3, 7, 11, 17, 31, 101]:
        t = modulus * h
        blocks, norms = {}, {}
        for label, (a, b) in FACTORS.items():
            c = constants[label]
            assert t % c == 0
            jj = 1 + a * (t // c), -b * (t // c)
            nn = norm(jj)
            assert mul(jj, (a, b)) == (t + a, b)
            assert nn * c == (t + a) ** 2 + b * b
            assert nn > 1 and primitive(jj)
            assert all(nn % p != 0 for p in ss)
            blocks[label], norms[label] = jj, nn
        assert all(gcd(x, y) == 1 for x, y in combinations(norms.values(), 2))
        rr = [gaussian_product([CORRECTIONS[i]] + [blocks[j] for j in INCIDENT[i]])
              for i in range(3)]
        assert [z[1] for z in rr] == [16, 90, 14]
        for i in range(3):
            assert primitive(rr[i]) and rr[i][0] > 0
            assert tuple(CONTENTS[i] * x for x in rr[i]) == polynomial_eval(rows[i], t)
            assert norm(rr[i]) == norm(CORRECTIONS[i]) * prod(norms[j] for j in INCIDENT[i])
        for (i, j), label, residual in [((0, 1), '12', -122),
                                       ((0, 2), '13', 18),
                                       ((1, 2), '23', 208)]:
            assert determinant(rr[i], rr[j]) == residual * norms['g'] * norms[label]
        b0 = gaussian_product(conjugate(z) for z in blocks.values())
        circle = [mul(circle_corrections[0], b0)]
        for i in range(3):
            bi = gaussian_product(z if label in INCIDENT[i] else conjugate(z)
                                  for label, z in blocks.items())
            circle.append(mul(circle_corrections[i + 1], bi))
        squared_radius = 876265 * prod(norms.values())
        assert all(norm(z) == squared_radius for z in circle)
        assert len(set(circle)) == 4
        circle_gcd = (0, 0)
        for point in circle:
            circle_gcd = gaussian_gcd(circle_gcd, point)
        assert norm(circle_gcd) == 1
        assert squared_radius * 4500 == prod((t + a) ** 2 + b * b for a, b in FACTORS.values())
        for i in range(3):
            assert mul(circle[i + 1], conjugate(rr[i])) == mul(circle[0], rr[i])
        denominator_lcm = (1, 0)
        for row in rr:
            denominator_lcm = gaussian_lcm(denominator_lcm, conjugate(row))
        assert norm(exact_divide(denominator_lcm, circle[0])) == 1
    print('Polynomial family, primitive progression, all determinants, and four-point circle verified.')


if __name__ == '__main__':
    main()
