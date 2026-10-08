#!/usr/bin/env python3
"""Exact rational contraction certificate for a full four-row real profile.

No numerical libraries are required.  The stored decimal roots and approximate
inverse are interpreted as exact rational numbers, not floating-point inputs.
"""

import json
from fractions import Fraction as Q
from pathlib import Path


def multiply(a, b):
    return (a[0] * b[0] - a[1] * b[1],
            a[0] * b[1] + a[1] * b[0])


def moments_and_jacobian(roots):
    powers = []
    for z in roots:
        row = [(Q(1), Q(0))]
        for _ in range(7):
            row.append(multiply(row[-1], z))
        powers.append(row)
    values, jacobian = [], []
    for i in range(4):
        masks = [j for j in range(15) if (j + 1) & (1 << i)]
        for k in range(1, 8):
            values.append(sum(powers[j][k][1] for j in masks))
            jacobian.append(
                [k * powers[j][k - 1][1] if j in masks else Q(0)
                 for j in range(14)]
                + [k * powers[j][k - 1][0] if j in masks else Q(0)
                   for j in range(14)])
    return values, jacobian


def infinity_norm(matrix):
    return max(sum(abs(x) for x in row) for row in matrix)


def main():
    source = Path(__file__).with_name("four_row_real_polynomial_certificate.json")
    data = json.loads(source.read_text())
    roots = [tuple(map(Q, row)) for row in data["roots"]]
    inverse = [list(map(Q, row)) for row in data["inverse"]]
    rho = Q(data["radius"])
    assert len(roots) == 15 and roots[-1] == (0, 1)
    assert len(inverse) == 28 and all(len(row) == 28 for row in inverse)
    assert rho == Q(1, 10**30)

    # The whole coordinate cube lies inside |z_T| < 2.
    assert all(abs(x) + abs(y) + 2 * rho < 2 for x, y in roots)
    values, jacobian = moments_and_jacobian(roots)
    af = [sum(a * f for a, f in zip(row, values)) for row in inverse]
    aj_error = [
        [Q(i == j) - sum(inverse[i][k] * jacobian[k][j]
                         for k in range(28))
         for j in range(28)]
        for i in range(28)]
    alpha = max(map(abs, af))
    beta = infinity_norm(aj_error)
    a_norm = infinity_norm(inverse)
    hessian_contribution = a_norm * 37632 * rho

    # These are checked with exact Fraction comparisons.
    assert alpha < Q(1, 10**60)
    assert beta < Q(1, 10**59)
    assert a_norm < 2700
    assert alpha <= rho / 2
    assert beta <= Q(1, 4)
    assert hessian_contribution <= Q(1, 4)
    contraction = beta + hessian_contribution
    assert contraction < Q(1, 2)
    assert alpha + contraction * rho < rho

    # Every root and conjugate is separated throughout the certified cube.
    doubled = roots + [(x, -y) for x, y in roots]
    separation = min(
        max(abs(a[0] - b[0]), abs(a[1] - b[1]))
        for i, a in enumerate(doubled) for b in doubled[i + 1:])
    assert separation - 2 * rho > Q(1, 50)
    assert min(abs(y) for _, y in roots) - rho > Q(1, 5)

    # The imaginary part of the constant term is the imaginary part of
    # the product of the eight roots.  Its movement is at most 1792*rho.
    constants = []
    for i in range(4):
        product = (Q(1), Q(0))
        for mask, z in enumerate(roots, 1):
            if mask & (1 << i):
                product = multiply(product, z)
        constants.append(product[1])
    assert min(constants) - 1792 * rho > Q(3, 1000)
    constant_separation = min(abs(a - b)
                              for i, a in enumerate(constants)
                              for b in constants[i + 1:])
    assert constant_separation - 3584 * rho > Q(1, 500)

    # Floats are used only for readable output after all exact checks pass.
    print("PASS: exact rational contraction and nondegeneracy certificate")
    print(f"||A F||_inf = {float(alpha):.6e}")
    print(f"||I-A J||_inf = {float(beta):.6e}")
    print(f"||A||_inf = {float(a_norm):.12f}")
    print(f"contraction bound = {float(contraction):.6e}")
    print(f"root/conjugate L-infinity separation > {float(separation-2*rho):.12f}")
    print("row imaginary constants:", [float(x) for x in constants])


if __name__ == "__main__":
    main()
