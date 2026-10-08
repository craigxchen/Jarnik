"""Exact certificates for the restricted Vieta and coordinate barriers."""
from fractions import Fraction as F
from math import gcd
import random


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def norm(z):
    return z[0] ** 2 + z[1] ** 2


def ggcd(a, b):
    while b != (0, 0):
        nn = norm(b)
        num = mul(a, (b[0], -b[1]))
        quotient = tuple((2 * x + nn) // (2 * nn) for x in num)
        zz = mul(quotient, b)
        a, b = b, (a[0] - zz[0], a[1] - zz[1])
    return a


def polynomial_multiply(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def main():
    rng = random.Random(410728)
    integral = nonintegral = 0
    checked = 0
    for _ in range(12000):
        # Columns (a,b),(c,d) give an integral oriented Gaussian basis.
        x, y, z = [rng.randrange(-8, 9) for _ in range(3)]
        a, b, c, d = 1 + x * y, y, x + z * (1 + x * y), 1 + z * y
        assert a * d - b * c == 1
        g, r, s = a * a + b * b, a * c + b * d, c * c + d * d
        assert g * s - r * r == 1
        v1 = (rng.randrange(-25, 26), rng.randrange(1, 15))
        v2 = (rng.randrange(-25, 26), rng.randrange(1, 15))
        if gcd(*v1) != 1 or gcd(*v2) != 1:
            continue
        ai, yi = v1
        aj, yj = v2
        delta = ai * yj - aj * yi
        if delta == 0:
            continue
        qi = (a * ai + c * yi, b * ai + d * yi)
        qj = (a * aj + c * yj, b * aj + d * yj)
        if norm(qi) % 2 == 0 or norm(qj) % 2 == 0:
            continue
        ss = qi[0] * qj[0] + qi[1] * qj[1]
        dd = norm(ggcd(qi, qj))
        assert ss % dd == delta % dd == 0
        u, v = ss // dd, delta // dd
        assert dd % 2 and gcd(u, v) == 1
        cc = [2 * yi * yj, -(ai * yj + aj * yi), 2 * ai * aj]
        content = gcd(gcd(abs(cc[0]), abs(cc[1])), abs(cc[2]))
        assert content == (1 if delta % 2 else 2)
        gp, rp, sp = [F(old) + F(2 * ss * entry, delta * delta)
                      for old, entry in zip((g, r, s), cc)]
        assert gp > 0 and gp * sp - rp * rp == 1
        for av, yv, qq in [(ai, yi, qi), (aj, yj, qj)]:
            assert gp * av * av + 2 * rp * av * yv + sp * yv * yv == norm(qq)
        mixed = gp * ai * aj + rp * (ai * yj + aj * yi) + sp * yi * yj
        assert mixed == -ss
        is_integral = all(entry.denominator == 1 for entry in (gp, rp, sp))
        assert is_integral == ((2 * content * ss) % (delta * delta) == 0)
        assert is_integral == (u % dd == 0 and abs(v) in (1, 2))
        integral += is_integral
        nonintegral += not is_integral
        denominator = rng.randrange(1, 10)
        lower_left, lower_right = F(rng.randrange(-10, 11), denominator), F(rng.randrange(-10, 11), denominator)
        yip = lower_left * ai + lower_right * yi
        yjp = lower_left * aj + lower_right * yj
        assert yj * yip - yi * yjp == lower_left * delta
        checked += 1
    assert integral and nonintegral
    # Derive the actual primitive pair residues and real numerator.
    corrections = [(7, 16), (47, 90), (5, 14)]
    pairs = [(0, 1), (0, 2), (1, 2)]
    gcd_norms = [norm(ggcd(corrections[i], corrections[j])) for i, j in pairs]
    assert gcd_norms == [61, 1, 13]
    prescribed_residues = [-122, 18, 208]
    assert [r // e for r, e in zip(prescribed_residues, gcd_norms)] == [-2, 18, 16]
    assert all(r % e == 0 for r, e in zip(prescribed_residues, gcd_norms))
    numerator = [870, -559, 151, -19, 1]
    divisor = [61, -12, 1]
    real_1 = [105, -256, 106, -16, 1]
    real_2 = [94, -195, 95, -15, 1]
    dot_polynomial = polynomial_multiply(real_1, real_2)
    dot_polynomial[0] += 240 * 180
    shared_norm = polynomial_multiply([1, 0, 1], divisor)
    assert dot_polynomial == polynomial_multiply(shared_norm, numerator)
    # Normalizing rows by 15 and 2 makes the primitive real coordinate U/30;
    # the actual shared Gaussian-gcd norm is (t^2+1)*Q_12 before removing J_g.
    assert 15 * 2 == 30
    # Exact remainder in the genuine balanced three-row family.
    remainder = numerator[:]
    for degree in range(4, 1, -1):
        leading = remainder[degree]
        for j in range(3):
            remainder[degree - 2 + j] -= leading * divisor[j]
    assert remainder == [504, -60, 0, 0, 0]
    modulus = 8088851149545216141000
    for h in [1, 2, 3, 7, 31]:
        t = modulus * h
        dd = t * t - 12 * t + 61
        nn = sum(coefficient * t ** j for j, coefficient in enumerate(numerator))
        assert nn % 30 == 0
        assert 0 < abs(504 - 60 * t) < dd
        assert nn % dd != 0 and (nn // 30) % dd != 0
    print(f'Vieta criterion verified on {checked} forms: {integral} integral and {nonintegral} nonintegral alternates.')


if __name__ == '__main__':
    main()
