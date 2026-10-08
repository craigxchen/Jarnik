"""Exact checks for the adaptive alignment height floor on the cubic family.

The note proves the inequalities for all admissible parameters; these finite
cases check the formulas and both independent least-radius constructions.
"""

from math import gcd

from check_endpoint_swap_content_and_overlap import (
    all_edge_radius,
    primitive_radius,
    reduce_denominator,
)
from check_least_radius_formula import (
    exact_div,
    gcd_gaussian,
    lcm_all,
    mul,
    norm,
)


K0 = 40 * 45360 * 60 * 25920


def family(t):
    P = t**4 + 25 * t**2 - 36
    Q = 60 * t
    B = (P, Q)
    U = mul(mul((t, 2), (t, 3)), (t, -6))
    assert mul(U, (t, 1)) == B
    d = 20 * (t**2 + 4) * (t**2 + 9) * (t**2 + 36)
    v = 40 * (t**2 + 18)
    assert d == 20 * norm(U)
    assert gcd(P, Q) == 36
    DB = reduce_denominator(B)
    assert norm(DB) == norm(B) // 36**2
    assert norm(DB) * 36**2 >= t**8
    other = (
        (48 * (t**4 + 10 * t**2 + 9), 12),
        (30 * (t**4 + 37 * t**2 + 36), 30),
        (d, v),
    )
    return B, other


def check_case(t, p, q, B, interiors):
    assert gcd(p, q) == 1 and p > 0 and q > 0 and p != q
    P, Q = B
    H = p + q
    d, v = interiors[-1]
    C = (q * d + p * v * P, p * v * Q)
    A = C[0]
    gamma = gcd(*C)
    assert A % t == (25920 * (q - p)) % t
    assert d % (t**2 + 18) == 45360
    assert gamma <= p * gcd(A, v) * gcd(A, Q)
    assert gcd(A, v) <= q * 40 * 45360
    assert gcd(A, Q) <= 60 * 25920 * abs(p - q)
    assert gamma <= K0 * p * q * abs(p - q)

    DC = reduce_denominator(C)
    assert norm(DC) * 2 * gamma**2 >= norm(C)
    assert norm(DC) * gamma**2 >= 200 * H**2 * t**12

    U = mul(mul((t, 2), (t, 3)), (t, -6))
    S = (t, 1)
    E = (
        20 * q * U[0] + p * v * t,
        -20 * q * U[1] + p * v,
    )
    assert mul(U, E) == C
    gSE = gcd_gaussian(S, E)
    exact_div((1200 * q, 0), gSE)
    assert norm(gcd_gaussian(B, C)) == norm(U) * norm(gSE)
    assert norm(gcd_gaussian(S, (1200 * q, 0))) == gcd(norm(S), 1200 * q)
    assert norm(gcd_gaussian(B, C)) <= 9600 * q * t**6

    subset = [(1, 0), B, C]
    target = [(1, 0), B] + [
        (q * di + p * vi * P, p * vi * Q) for di, vi in interiors
    ]
    for (di, vi), row in zip(interiors, target[2:]):
        g = gcd(di, vi)
        d0, v0 = di // g, vi // g
        a, b = gcd(p, d0), gcd(q, v0)
        alpha, beta = (q // b) * (d0 // a), (p // a) * (v0 // b)
        assert gcd(alpha, beta) == 1
        assert gcd(*row) == g * a * b * gcd(alpha + beta * P, Q)
    for i in range(3):
        for j in range(i + 1, 3):
            di, vi = interiors[i]
            dj, vj = interiors[j]
            W = vj * di - vi * dj
            overlap = gcd_gaussian(target[2 + i], target[2 + j])
            exact_div(mul((W, 0), gcd_gaussian(B, (q, 0))), overlap)
    subset_lcm = norm(lcm_all([reduce_denominator(B), DC]))
    subset_n = primitive_radius(subset)
    full_n = primitive_radius(target)
    assert subset_n == subset_lcm == all_edge_radius(subset)
    assert full_n == all_edge_radius(target)
    assert full_n % subset_n == 0

    assert full_n * 62208 * gamma**2 * q >= t**14 * H**2
    assert (full_n * 62208 * K0**2 * p**2 * q**3 * (p - q)**2
            >= t**14 * H**2)
    assert full_n * 62208 * K0**2 * H**3 * (p - q)**2 >= t**14
    assert full_n * 62208 * K0**2 * H**5 >= t**14


def check_shared_endpoint_prime():
    t, p, q = 600, 170, 169
    B, _ = family(t)
    D1 = 4 * (t**2 + 1) * (t**2 + 9)
    D2 = (t**2 + 1) * (t**2 + 36)
    rows = [B] + [(q * D + p * B[0], p * B[1]) for D in (D1, D2)]
    pi, conjugate_pi = (3, -2), (3, 2)
    for row in rows:
        assert norm(gcd_gaussian(row, mul(pi, pi))) == 13
        assert norm(gcd_gaussian(row, conjugate_pi)) == 1
    denominators = [reduce_denominator(row) for row in rows]
    assert norm(gcd_gaussian(lcm_all(denominators),
                             mul(conjugate_pi, conjugate_pi))) == 13
    assert D1 - D2 == 3 * t**2 * (t**2 + 1)
    assert (D1 - D2) % 13 != 0


def main():
    count = 0
    for t in (600, 1200, 2400):
        B, interiors = family(t)
        for p, q in ((2, 1), (1, 2), (3, 2), (2, 3),
                     (t + 1, t), (t - 1, t), (t // 10 + 1, 1)):
            if gcd(p, q) != 1:
                continue
            check_case(t, p, q, B, interiors)
            count += 1
    assert count >= 18
    check_shared_endpoint_prime()
    print(f"adaptive height: {count} exact five-row cases and one shared-prime fixture passed")


if __name__ == "__main__":
    main()
