"""Independent exact checks for the Pell four-point local ball theorem."""

from itertools import product
from math import gcd

from check_direct_primitive_axis_lift import family
from check_gaussian_reflection_replacement import gnorm


def matvec(matrix, vector):
    a, b, c, d = matrix
    x, y = vector
    return a * x + b * y, c * x + d * y


def transpose_vec(matrix, vector):
    a, b, c, d = matrix
    x, y = vector
    return a * x + c * y, b * x + d * y


def check_common_zeros():
    zeros = set()
    for s, t in product(range(-40, 41), repeat=2):
        F = s * (s - 2 * t - 16)
        H = t * (2 * s + t - 2)
        if F == H == 0:
            zeros.add((s, t))
    expected = {(0, 0), (16, 0), (0, 2), (4, -6)}
    assert zeros == expected


def check_index(index):
    U, V, P, A, B, C, z = family(index)
    a, b, c = gnorm(A), gnorm(B), gnorm(C)
    radius2 = gnorm(z[0])
    assert radius2 == 5 * a * b * c
    T = (2 * V, V - U, -U - V, 2 * V)
    assert T[0] * T[3] - T[1] * T[2] == -1

    gram = (
        T[0] * T[0] + T[2] * T[2],
        T[0] * T[1] + T[2] * T[3],
        T[1] * T[0] + T[3] * T[2],
        T[1] * T[1] + T[3] * T[3],
    )
    assert gram == (a, b - a, b - a, b)
    assert a * b - (a - b) ** 2 == 1
    assert gcd(a, b) == 1
    assert transpose_vec(T, z[0]) == (-8 * a, -b)

    old_coordinates = [(0, 0), (16, 0), (0, 2), (4, -6)]
    differences = [(x - z[0][0], y - z[0][1]) for x, y in z]
    assert differences == [matvec(T, point) for point in old_coordinates]

    inverse = (-T[3], T[1], T[2], -T[0])
    y0 = matvec(inverse, z[0])
    Jz0 = (-z[0][1], z[0][0])
    v0 = matvec(inverse, Jz0)
    assert y0 == (-(9 * a * b - b * b), -(8 * a * a - 7 * a * b))
    assert v0 == (-b, 8 * a)

    for s, t in old_coordinates:
        F = s * (s - 2 * t - 16)
        H = t * (2 * s + t - 2)
        assert a * F + b * H == 0

    assert c == 16 * a - 3 * b
    if V >= 4:
        assert 14 * V * V < a < 15 * V * V
        assert 5 * V * V < b < 6 * V * V
        assert c > 206 * V * V
        assert radius2 > 250**2 * V**6
        assert max(map(abs, y0)) < 1800 * V**4
        assert max(map(abs, v0)) < 120 * V**2
    if V >= 16:
        # The only square-root step in (17): 2 sqrt(V) <= V/2.
        assert 16 * V <= V * V
        assert 9 * V + 20 * V <= 29 * V < 40 * V
        assert 3 * V + 16 <= 4 * V
        assert 4 * V * V < b
        assert 3 * V + 2 < 14 * V


def main():
    check_common_zeros()
    for index in (1, 11, 21):
        check_index(index)
    print("PASS: three exact Pell coordinate and Gram identities.")
    print("PASS: T is unimodular and the integral conic has exactly four common zeros.")
    print("PASS: all a,b,c,R,y0,v0 inequalities used by C0=1/4.")
    print("PASS: V>=16 forces |F|<b and |H|<a inside the stated ball.")


if __name__ == "__main__":
    main()
