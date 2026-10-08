"""Exact local fixtures for cross-support Hankel conditions; not eight-row solutions."""

from fractions import Fraction as F


def solve(matrix, rhs):
    a = [list(map(F, row)) + [F(value)]
         for row, value in zip(matrix, rhs)]
    for i in range(len(a)):
        pivot = next(j for j in range(i, len(a)) if a[j][i])
        a[i], a[pivot] = a[pivot], a[i]
        scale = a[i][i]
        a[i] = [value / scale for value in a[i]]
        for j in range(len(a)):
            if j != i:
                scale = a[j][i]
                a[j] = [u - scale * v for u, v in zip(a[j], a[i])]
    return [row[-1] for row in a]


def polynomial(roots):
    coeff = [F(1)]
    for root in roots:
        updated = [F(0)] * (len(coeff) + 1)
        for i, value in enumerate(coeff):
            updated[i] -= root * value
            updated[i + 1] += value
        coeff = updated
    return coeff


def power_sums(coeff, count):
    """Newton sums for a monic polynomial with coefficients in ascending order."""
    degree = len(coeff) - 1
    assert coeff[-1] == 1
    out = [F(degree)]
    for k in range(1, count + 1):
        value = sum(coeff[degree - j] * out[k - j]
                    for j in range(1, min(k - 1, degree) + 1))
        if k <= degree:
            value += k * coeff[degree - k]
        out.append(-value)
    return out


def check(roots, expected_first):
    d = len(roots)
    s = (d - 1) // 2
    x = [F(z * z) for z in roots]
    v = [F(1, z) for z in roots]
    assert len(set(x)) == d and min(x) > 0
    assert all(sum(z ** (2 * r + 1) for z in roots) == 0
               for r in range(s))
    coeff = polynomial(roots)
    assert coeff[0] != 0
    assert all(coeff[i] == 0 for i in range(2, d, 2))
    assert sum(v) ** 2 == sum(t * t for t in v)
    h_coeff = coeff[1::2]
    squared_h = [F(0)] * (2 * s + 1)
    for i, left in enumerate(h_coeff):
        for j, right in enumerate(h_coeff):
            squared_h[i + j] += left * right
    assert polynomial(x) == [-coeff[0] ** 2] + squared_h
    xi_powers = power_sums(h_coeff, 2 * s)
    moments = [sum(t ** k for t in x) for k in range(2 * s + 1)]
    assert all(moments[k] == 2 * xi_powers[k]
               for k in range(1, 2 * s + 1))
    for r1 in range(1, 2 * s):
        for r2 in range(1, 2 * s - r1 + 1):
            gap = moments[r1] * moments[r2] - 2 * moments[r1 + r2]
            assert gap > 0 if s >= 2 else gap == 0
    if s >= 2:
        largest, second = sorted(x, reverse=True)[:2]
        assert (largest ** s + (d - 1) * second ** s) ** 2 > (
            2 * largest ** (2 * s))
    values = []
    for k in range(1, s + 1):
        moment = lambda r: sum(t ** r for t in x)
        h = [moment(i) for i in range(1, k + 1)]
        matrix = [[moment(i + j) for j in range(1, k + 1)]
                  for i in range(1, k + 1)]
        projection = solve(matrix, h)
        fitted_norm = sum(a * b for a, b in zip(h, projection))
        residual = [1 - sum(projection[j - 1] * t ** j
                            for j in range(1, k + 1)) for t in x]
        assert all(sum(p * t ** j for p, t in zip(residual, x)) == 0
                   for j in range(1, k + 1))
        assert sum(p * p for p in residual) == d - fitted_norm
        if k == s:
            assert fitted_norm == d - 1
            assert all(w == sum(v) * p for w, p in zip(v, residual))
        else:
            assert fitted_norm < d - 1
        values.append(fitted_norm)
    assert values[0] == expected_first
    print(f"PASS: d={d}, fitted norms {values}; square identity and quadrature")


if __name__ == "__main__":
    check([-3, 1, 2], F(2))
    check([1, 5, 9, -7, -8], F(1100, 311))
