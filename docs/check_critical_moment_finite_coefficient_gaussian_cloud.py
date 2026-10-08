"""Exact rational fixture for the finite-coefficient Gaussian cloud."""

from fractions import Fraction


def multiply(a, b):
    c = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return c


def evaluate(coefficients, x):
    value = Fraction(0)
    for coefficient in reversed(coefficients):
        value = value * x + coefficient
    return value


def main():
    # E_0(t)=1+t-t^3/27=1-2T_3(t/6).
    base = [Fraction(1), Fraction(1), Fraction(0), -Fraction(1, 27)]
    q = Fraction(1, 10000)
    cloud = [Fraction(1), Fraction(0), -q]
    # Truncate E_0/(1-q t^2) at degree three.
    macro = [Fraction(1), Fraction(1), q, q - Fraction(1, 27)]
    full = multiply(cloud, macro)
    assert full[:4] == base
    assert full[2] == 0
    # For E=prod(1-y_j t), p_1=-[t]E and p_2=p_1^2-2[t^2]E.
    p1 = -full[1]
    p2 = p1 * p1 - 2 * full[2]
    assert p1 == -1 and p2 == 1
    assert evaluate(macro, -5) > 0 and evaluate(macro, -4) < 0
    assert evaluate(macro, -2) < 0 and evaluate(macro, -1) > 0
    assert evaluate(macro, 5) > 0 and evaluate(macro, 6) < 0
    # Three simple nonzero real macro roots have |t|>1; cloud roots
    # are exactly +/-100, with reciprocal-square mass 2/10000.
    assert evaluate(cloud, -100) == evaluate(cloud, 100) == 0
    assert 2 * q == Fraction(1, 5000)
    print("PASS: degree-five rational real-rooted product, e2=0, p1=-1, p2=1, cloud mass=1/5000")


if __name__ == "__main__":
    main()
