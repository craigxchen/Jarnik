"""Exact checks for the ordered terminal-circuit sign countermodel."""

from collections import defaultdict
from itertools import combinations


def linear_circuit(triple, slopes):
    a, b, c = triple
    return {a: slopes[c] - slopes[b],
            b: -(slopes[c] - slopes[a]),
            c: slopes[b] - slopes[a]}


def multiply_linear(left, right, scale=1):
    out = defaultdict(int)
    for i, a in left.items():
        for j, b in right.items():
            out[tuple(sorted((i, j)))] += scale * a * b
    return dict(out)


def determinant3(rows):
    a, b, c = rows
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def triangle(rows, triple):
    return determinant3([(rows[i][0] ** 2,
                          rows[i][0] * rows[i][1],
                          rows[i][1] ** 2) for i in triple])


def check(N):
    slopes = tuple(range(6))
    f = multiply_linear(linear_circuit((0, 1, 2), slopes),
                        linear_circuit((3, 4, 5), slopes), 9)
    second = multiply_linear(linear_circuit((0, 1, 3), slopes),
                             linear_circuit((2, 4, 5), slopes), -1)
    for monomial, coefficient in second.items():
        f[monomial] = f.get(monomial, 0) + coefficient
    f = {monomial: coefficient for monomial, coefficient in f.items()
         if coefficient}

    expected = {
        (0, 2): -2, (0, 3): 9, (0, 4): -12, (0, 5): 5,
        (1, 2): 3, (1, 3): -18, (1, 4): 27, (1, 5): -12,
        (2, 3): 8, (2, 4): -18, (2, 5): 9, (3, 4): 3, (3, 5): -2,
    }
    assert f == expected
    assert any(c > 0 for c in f.values())
    assert any(c < 0 for c in f.values())

    # D_1 f = D_b f = 0, equivalently every terminal circuit vanishes.
    for i in range(6):
        assert sum(f.get(tuple(sorted((i, j))), 0)
                   for j in range(6) if j != i) == 0
        assert sum(slopes[j] * f.get(tuple(sorted((i, j))), 0)
                   for j in range(6) if j != i) == 0

    # Reconstruction at the actual point x_i=b_i^2.
    reconstructed = sum(c * slopes[i] ** 2 * slopes[j] ** 2
                        for (i, j), c in f.items())
    assert reconstructed == 0

    # For binary rows (N,i), Q=9*T012*T345-T013*T245 is zero.
    rows = tuple((N, i) for i in range(6))
    value = (9 * triangle(rows, (0, 1, 2)) * triangle(rows, (3, 4, 5))
             - triangle(rows, (0, 1, 3)) * triangle(rows, (2, 4, 5)))
    assert value == 0


for N in (1, 7, 101):
    check(N)

print("PASS: ordered circuits have mixed signs, vanish on reconstruction, and arise from a nonzero terminal coherent array.")
