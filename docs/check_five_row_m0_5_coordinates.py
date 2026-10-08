"""Exact M_0,5 chart from O,P,3P,7P on y^2=x^3-2x."""

from fractions import Fraction as F


O = None
P = (F(2), F(2))


def add(Q, R):
    if Q is None:
        return R
    if R is None:
        return Q
    x1, y1 = Q
    x2, y2 = R
    if x1 == x2 and y1 == -y2:
        return None
    if Q == R:
        slope = (3 * x1 * x1 - 2) / (2 * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope * slope - x1 - x2
    y3 = -(y1 + slope * (x3 - x1))
    return x3, y3


def multiple(n, Q=P):
    if n < 0:
        return multiple(-n, (Q[0], -Q[1]))
    out = None
    base = Q
    while n:
        if n & 1:
            out = add(out, base)
        base = add(base, base)
        n >>= 1
    return out


def line_slope_from_P(Q):
    if Q is None:
        return None
    if Q == P:
        return F(5, 2)  # tangent at P
    x, y = Q
    if x == P[0]:
        return None  # vertical line through P and -P
    return (y - P[1]) / (x - P[0])


C = multiple(3)
D = multiple(7)
xP = P[0]
xC = C[0]
xD = D[0]
sC = line_slope_from_P(C)
sD = line_slope_from_P(D)


def check_exact_constants():
    assert C == (F(338), F(6214))
    assert xD == F(9374597220578, 5627138321281)
    assert sC == F(1553, 84)
    assert sD == F(927358740467713, 358778440454952)


def a_coordinate(Q):
    """Pencil through O, normalized AC -> 0, AD -> 1, AB -> infinity."""
    if Q is None:
        # The limiting value at O is the tangent-line direction.
        return (xD - xP) / (xD - xC)
    x, _ = Q
    if x == xP:
        return None
    return (x - xC) * (xD - xP) / ((x - xP) * (xD - xC))


def b_coordinate(Q):
    """Pencil through P, normalized BC -> 0, BD -> 1, BA -> infinity."""
    slope = line_slope_from_P(Q)
    if slope is None:
        return None
    return (slope - sC) / (sD - sC)


def check_normalization():
    assert a_coordinate(C) == 0
    assert a_coordinate(D) == 1
    assert a_coordinate(P) is None
    assert b_coordinate(C) == 0
    assert b_coordinate(D) == 1
    assert b_coordinate(O) is None
    assert a_coordinate(multiple(-1)) is None
    assert b_coordinate(multiple(-1)) is None


def check_boundary_indices():
    # Exceptional points O,P,3P,7P and third intersections -(a+b)P of
    # the six lines through pairs of centers.
    indices = {0, 1, 3, 7, -1, -3, -7, -4, -8, -10}
    assert len(indices) == 10
    assert len(indices) == len({multiple(k) for k in indices})


def check_boundary_evaluations():
    # The labels in the chart are x_1=infinity, x_2=0, x_3=1,
    # x_4=a, x_5=b.  A value at infinity is represented by None.
    expected = {
        0: (("finite",), ("inf",)),
        1: (("inf",), ("finite",)),
        3: ((0,), (0,)),
        7: ((1,), (1,)),
        -1: (("inf",), ("inf",)),
        -3: ((0,), ("finite",)),
        -7: ((1,), ("finite",)),
        -4: (("finite",), (0,)),
        -8: (("finite",), (1,)),
        -10: (("equal",), ("equal",)),
    }

    def category(value):
        if value is None:
            return "inf"
        return value

    for k, (a_expected, b_expected) in expected.items():
        a, b = a_coordinate(multiple(k)), b_coordinate(multiple(k))
        if a_expected == ("finite",):
            assert a is not None and a not in (0, 1)
        elif a_expected == ("equal",):
            assert a == b
        else:
            assert category(a) == a_expected[0]
        if b_expected == ("finite",):
            assert b is not None and b not in (0, 1)
        elif b_expected == ("equal",):
            assert a == b
        else:
            assert category(b) == b_expected[0]


def check_chart_away_from_boundary():
    boundary = {0, 1, 3, 7, -1, -3, -7, -4, -8, -10}
    samples = {}
    for n in (2, 4, 5, 6, 8, 9, 11):
        assert n not in boundary
        Q = multiple(n)
        a, b = a_coordinate(Q), b_coordinate(Q)
        assert a is not None and b is not None
        assert a not in (0, 1) and b not in (0, 1) and a != b
        samples[n] = (a, b)
    return samples


def check_boundary_chart_dictionary():
    # In the chart (infinity,0,1,a,b), the seven visible boundaries are
    # a=0:D24, b=0:D25, a=1:D34, b=1:D35, a=infty:D14,
    # b=infty:D15, and a=b:D45.  The three blown-up centers are
    # (0,0):D13, (1,1):D12, (infty,infty):D23.
    dictionary = {
        "a=0": "D24", "b=0": "D25", "a=1": "D34", "b=1": "D35",
        "a=inf": "D14", "b=inf": "D15", "a=b": "D45",
        "(0,0)": "D13", "(1,1)": "D12", "(inf,inf)": "D23",
    }
    assert len(dictionary) == 10
    assert set(dictionary.values()) == {
        "D12", "D13", "D14", "D15", "D23", "D24", "D25", "D34", "D35", "D45"
    }


def main():
    check_exact_constants()
    check_normalization()
    check_boundary_indices()
    check_boundary_evaluations()
    samples = check_chart_away_from_boundary()
    check_boundary_chart_dictionary()
    print("Exact pencil normalization and ten boundary indices pass.")
    print("a(Q)=((x-x3)(x7-xP))/((x-xP)(x7-x3)); b(Q)=((s-s3)/(s7-s3)), s=(y-yP)/(x-xP).")
    print("Boundary dictionary: D24:a=0, D25:b=0, D34:a=1, D35:b=1, D14:a=inf, D15:b=inf, D45:a=b.")
    print("Blown-up centers: (0,0)=D13, (1,1)=D12, (inf,inf)=D23.")
    for n, (a, b) in samples.items():
        print(f"N={n}: a={a}, b={b}")


if __name__ == "__main__":
    main()
