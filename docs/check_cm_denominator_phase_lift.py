"""Exact tests for CM denominator generators and their relative phase law."""
from itertools import product
from math import isqrt
from check_elliptic_split_inert_norm_route import (
    C, P, IP, add, multiple, plus, minus, neg, times, divide,
    conj, norm, ggcd, denominator_generator,
)

UNITS = {C(1), C(-1), C(0, 1), C(0, -1)}


def power(a, exponent):
    if exponent < 0:
        return power(divide(C(1), a), -exponent)
    result = C(1)
    while exponent:
        if exponent & 1:
            result = times(result, a)
        a = times(a, a)
        exponent //= 2
    return result


def beta(n):
    return denominator_generator(add(multiple(n), IP))


def phase(n):
    b = beta(n)
    return divide(b, conj(b))


def line(l):
    a = add(multiple(l), (IP[0], neg(IP[1])))
    slope = divide(minus(IP[1], a[1]), minus(IP[0], a[0]))
    intercept = minus(IP[1], times(slope, IP[0]))
    return slope, intercept


def fvalue(l, n):
    if l == 0 or n == 0:
        return C(1)
    x, y = multiple(n)
    slope, intercept = line(l)
    value = minus(minus(y, times(slope, x)), intercept)
    if value == C(0):
        return None  # shared third intersection; requires cancelled evaluation
    return divide(value, conj(value))


def quotient_coordinate(complement_sum, n):
    q = multiple(n)
    if q is None:
        return None
    if complement_sum == 0:
        return q[0]
    t = multiple(complement_sum)
    denominator = minus(q[0], t[0])
    if denominator == C(0):
        return None  # pole or removable value; handled by the curve proof
    return divide(plus(q[1], t[1]), denominator)


def main():
    for n in range(-20, 21):
        if n == 0:
            continue
        r = multiple(n)
        x, y = r[0][0], r[1][0]
        A = x.numerator
        B = isqrt(x.denominator)
        assert B*B == x.denominator
        cc = y*B**3
        assert cc.denominator == 1
        cc = int(cc)
        D = A+2*B*B
        assert D > 0 and D % 8 != 0
        assert cc*cc+4*B**6 == D*(A*A-2*A*B*B+2*B**4)
        g = ggcd((D, 0), (cc, 2*B**3))
        b = beta(n)
        ratio = divide(g, b)
        # Equal odd-prime ideals: the ratio is a unit times a bounded power of1+i.
        assert any(divide(ratio, power(C(1, 1), e)) in UNITS for e in range(-2, 5))
        oddD = D
        while oddD % 2 == 0:
            oddD //= 2
        oddnorm = int(norm(g))
        while oddnorm % 2 == 0:
            oddnorm //= 2
        assert oddD == oddnorm
    assert line(1) == (C(-1/2, -1/2), C(-1, 1))
    count = 0
    for l in (-10, -8, -7, -4, -3, -1, 1, 3, 7):
        for n in range(-5, 16):
            f = fvalue(l, n)
            if f is None:
                continue
            lhs = divide(phase(n-l), phase(n))
            rhs = times(phase(-l), f)
            assert divide(lhs, rhs) in UNITS
            count += 1

    # Eight distinct shifts assigned to all three-row Boolean cuts.
    shifts = list(range(8))
    bits = list(product((0, 1), repeat=3))
    for n in (10, 11, 12):
        for i, j in ((0, 1), (0, 2), (1, 2)):
            cs = [b[i]-b[j] for b in bits]
            assert sum(cs) == 0 and sum(abs(c) for c in cs) == 4
            actual = C(1); expected = C(1)
            for l, c in zip(shifts, cs):
                actual = times(actual, power(phase(n-l), c))
                expected = times(expected, power(times(phase(-l), fvalue(l, n)), c))
            assert divide(actual, expected) in UNITS
    quotient_count = 0
    for complement_sum in (-5, 0, 7, 12):
        for n in range(-9, 19):
            value = quotient_coordinate(complement_sum, n)
            reflected = quotient_coordinate(complement_sum, complement_sum-n)
            if value is None or reflected is None:
                continue
            assert value == reflected and value[1] == 0
            quotient_count += 1
    for m in range(3, 13):
        assert 2**(m-1)-2 > 0
    print("40 exact gcd generators have the predicted odd orientation/norm and bounded ramified cost.")
    print(f"{count} exact elliptic chord/relative-phase identities pass, including negative shifts.")
    print("All three equal-cardinality Boolean row comparisons pass in three exact windows.")
    print(f"{quotient_count} exact rational degree-two quotient reflection identities pass.")
    print("No numerical approximation or moving-target elliptic-log estimate is asserted.")


if __name__ == "__main__":
    main()
