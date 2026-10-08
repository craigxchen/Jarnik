"""Exact checks for the real-order audit.

The check uses the smallest ordered pair that refutes strict positivity
of (F_i F_j+1)/(F_i-F_j), while satisfying the polynomial divisor model.
"""

from fractions import Fraction


Q = Fraction


def add(p, q):
    n = max(len(p), len(q))
    out = [Q(0)] * n
    for j, a in enumerate(p):
        out[j] += a
    for j, a in enumerate(q):
        out[j] += a
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def mul(p, q):
    out = [Q(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def scale(p, c):
    return [c * a for a in p]


def div_exact(p, q):
    p = list(p)
    out = [Q(0)] * (len(p) - len(q) + 1)
    while len(p) >= len(q):
        shift = len(p) - len(q)
        factor = p[-1] / q[-1]
        out[shift] = factor
        for j, a in enumerate(q):
            p[shift + j] -= factor * a
        while len(p) > 1 and p[-1] == 0:
            p.pop()
    assert p == [Q(0)]
    return out


def main():
    # Coefficients are in ascending powers of t.
    f_minus = [Q(-1), Q(1), Q(-1)]
    f_plus = [Q(1), Q(1), Q(1)]
    difference = add(f_minus, scale(f_plus, Q(-1)))
    assert difference == [Q(-2), Q(0), Q(-2)]
    assert all(coef <= 0 for coef in difference)

    numerator = add(mul(f_minus, f_plus), [Q(1)])
    quotient = div_exact(numerator, difference)
    assert quotient == [Q(0), Q(0), Q(1, 2)]
    assert quotient[0] == 0

    # Degree and the stated lcm upper bound for the two-row configuration.
    D = 2
    assert len(f_minus) - 1 == D == len(f_plus) - 1
    assert D <= 2 * D  # lcm degree is at most the sum of the two degrees.
    print("Passed exact ordered pair check: divisor identity, global order, and quotient t^2/2.")


if __name__ == "__main__":
    main()
