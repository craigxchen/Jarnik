"""Exact odd-multiple configurations on the disjoint elliptic six-curve.

The quotient quartic is used only as an exact computational model.  The
three inverse discriminants being squares certifies that the resulting
target value lifts to the rational multiquadratic source.
"""
from fractions import Fraction as F
from math import gcd, isqrt


A = 67144321
B = 212248804
C = 46513
D = 208624
K = -4 * B * D
LAMBDAS = (F(0), F(-240, 5329), F(8, 257))
ANCHORS = ((1, 0), (0, 1), (1, 1))  # infinity, 0, 1


def cubic_constants():
    roots = (4 * A * D, B * D, 4 * B * C)
    return (sum(roots),
            roots[0] * roots[1] + roots[0] * roots[2] + roots[1] * roots[2],
            roots[0] * roots[1] * roots[2])


def cubic_add(left, right):
    a2, a4, _ = cubic_constants()
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2 and y1 == -y2:
        return None
    if left == right:
        slope = (3 * x1 * x1 + 2 * a2 * x1 + a4) / (2 * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope * slope - a2 - x1 - x2
    y3 = -y1 + slope * (x1 - x3)
    return x3, y3


def cubic_neg(point):
    return None if point is None else (point[0], -point[1])


def cubic_multiple(n, point):
    result = None
    addend = point
    while n:
        if n & 1:
            result = cubic_add(result, addend)
        addend = cubic_add(addend, addend)
        n //= 2
    return result


def quartic_value(t):
    return t * (A * t - B) * (t - 4) * (C * t - D)


def rational_sqrt(value):
    assert value >= 0
    numerator = isqrt(value.numerator)
    denominator = isqrt(value.denominator)
    assert numerator * numerator == value.numerator
    assert denominator * denominator == value.denominator
    return F(numerator, denominator)


def inverse_x_roots(lam, target):
    # This is N_lambda(x)-T D_lambda(x)=0 before the X normalization.
    aa = 1 - lam * target
    bb = -73 * lam - target + 11 * lam * target
    cc = -60 * lam + target + 38 * lam * target
    square_root = rational_sqrt(bb * bb - 4 * aa * cc)
    return ((-bb + square_root) / (2 * aa),
            (-bb - square_root) / (2 * aa))


def normalized_x(x):
    return (x + 1) / (16 - 4 * x)


def primitive_row(value):
    numerator, denominator = value.numerator, value.denominator
    common = gcd(abs(numerator), denominator)
    return numerator // common, denominator // common


def config_for_odd_multiple(multiple):
    # T=K/x and Y=y*K/x^2 under the stated quartic/cubic conversion.
    r0_t = F(-1, 2)
    r0_y = F(45299439, 4)
    r0 = (K / r0_t, K * r0_y / r0_t**2)
    cubic_point = cubic_multiple(multiple, r0)
    target = K / cubic_point[0]
    quartic_y = cubic_point[1] * K / cubic_point[0]**2
    assert quartic_y * quartic_y == quartic_value(target)

    moving = []
    for lam in LAMBDAS:
        old_roots = inverse_x_roots(lam, target)
        assert all(
            (x * x - 73 * lam * x - 60 * lam)
            / (lam * x * x + (1 - 11 * lam) * x - 1 - 38 * lam)
            == target
            for x in old_roots
        )
        moving.append(tuple(normalized_x(x) for x in old_roots))
    rows = list(ANCHORS)
    rows.extend(primitive_row(value) for pair in moving for value in pair)
    assert len(set(rows)) == 9
    assert all(gcd(abs(row[0]), row[1]) == 1 for row in rows)
    return target, quartic_y, moving, rows


def check():
    # q=[2]P-R0 for P=nR0, so n=2,3 give q=3R0,5R0.
    r0_t = F(-1, 2)
    r0_y = F(45299439, 4)
    r0 = (K / r0_t, K * r0_y / r0_t**2)
    a2, a4, a6 = cubic_constants()
    assert r0[1] * r0[1] == r0[0]**3 + a2 * r0[0]**2 \
        + a4 * r0[0] + a6
    for n, odd_multiple in ((2, 3), (3, 5)):
        point = cubic_multiple(n, r0)
        assert cubic_add(cubic_add(point, point), cubic_neg(r0)) \
            == cubic_multiple(odd_multiple, r0)
        target, quartic_y, moving, rows = config_for_odd_multiple(odd_multiple)
        assert quartic_y
        selected = list(ANCHORS) + [primitive_row(pair[0]) for pair in moving]
        assert len(set(selected)) == 6
        print(f"n={n}, q={odd_multiple}R0")
        print(f"  T={target}")
        print(f"  Y={quartic_y}")
        for lam, pair in zip(LAMBDAS, moving):
            print(f"  lambda={lam}: X={pair}")
        print(f"  selected six-point rows={selected}")
        print(f"  primitive rows={rows}")
        # The three square roots are an independent exact lift check.
        checks = (
            400 * target * (target - 4),
            F(400, 73**4) * target * (67144321 * target - 212248804),
            F(400, 257**2) * (target - 4)
            * (46513 * target - 208624),
        )
        assert all(rational_sqrt(value) for value in checks)
        # These three anchor determinants are the expected primitive units.
        assert rows[0][0] * rows[1][1] - rows[1][0] * rows[0][1] == 1
        assert rows[1][0] * rows[2][1] - rows[2][0] * rows[1][1] == -1
        assert rows[2][0] * rows[0][1] - rows[0][0] * rows[2][1] == -1


if __name__ == "__main__":
    check()
    print("Exact odd-multiple source lifts and primitive directions pass.")
