#!/usr/bin/env python3
"""Exact finite checks for composite_inert_crt_collision_budget.md."""

from __future__ import annotations

from itertools import combinations, product
from math import atan, comb, gcd, prod, sqrt


Gaussian = tuple[int, int]


def gmul(a: Gaussian, b: Gaussian) -> Gaussian:
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def gnorm(a: Gaussian) -> int:
    return a[0]*a[0]+a[1]*a[1]


def grem(a: Gaussian, b: Gaussian) -> Gaussian:
    den = gnorm(b)
    real = a[0]*b[0]+a[1]*b[1]
    imag = a[1]*b[0]-a[0]*b[1]
    qr = (2*real+den)//(2*den)
    qi = (2*imag+den)//(2*den)
    return (a[0]-qr*b[0]+qi*b[1], a[1]-qr*b[1]-qi*b[0])


def ggcd(a: Gaussian, b: Gaussian) -> Gaussian:
    while b != (0, 0):
        a, b = b, grem(a, b)
    return a


def gquot(a: Gaussian, b: Gaussian) -> Gaussian:
    den = gnorm(b)
    real = a[0]*b[0]+a[1]*b[1]
    imag = a[1]*b[0]-a[0]*b[1]
    assert real % den == 0 and imag % den == 0
    return (real//den, imag//den)


def e(m: int, q: int) -> int:
    b, r = divmod(m, q)
    return q*comb(b, 2)+r*b


def factor(n: int) -> dict[int, int]:
    answer = {}
    divisor = 2
    while divisor*divisor <= n:
        while n % divisor == 0:
            answer[divisor] = answer.get(divisor, 0)+1
            n //= divisor
        divisor += 1
    if n > 1:
        answer[n] = answer.get(n, 0)+1
    return answer


def conic_class_count(modulus: int, norm: int) -> int:
    return sum(1 for x in range(modulus) for y in range(modulus)
               if (x*x+y*y-norm) % modulus == 0)


def build_points(m: int, n: int) -> tuple[list[Gaussian], int]:
    hs = [(2*(n+j), 1) for j in range(m)]
    points = []
    for j in range(m):
        z = hs[j]
        for ell in range(m):
            if ell != j:
                z = gmul(z, (hs[ell][0], -hs[ell][1]))
        points.append(z)
    norms = {gnorm(z) for z in points}
    assert len(norms) == 1
    return points, norms.pop()


def cofactor(a: Gaussian, b: Gaussian) -> Gaussian:
    difference = (a[0]-b[0], a[1]-b[1])
    return gquot(difference, ggcd(a, b))


def check_actual_families() -> None:
    inert_moduli = [3, 7, 9, 11, 21, 27, 33, 49, 63, 77, 99, 147]
    comparisons = 0
    for m in (8, 11, 15, 22):
        for n in (1, 4, 17):
            points, common_norm = build_points(m, n)
            cofactors = {(i, j): cofactor(points[i], points[j])
                         for i, j in combinations(range(m), 2)}
            for modulus in inert_moduli:
                assert gcd(modulus, common_norm) == 1
                collisions = 0
                for i, j in combinations(range(m), 2):
                    direct = ((points[i][0]-points[j][0]) % modulus == 0
                              and (points[i][1]-points[j][1]) % modulus == 0)
                    primitive = (cofactors[i, j][0] % modulus == 0
                                 and cofactors[i, j][1] % modulus == 0)
                    expected = (i-j) % modulus == 0
                    assert direct == primitive == expected
                    collisions += expected
                    comparisons += 1
                assert collisions == e(m, modulus)
    print(f"checked {comparisons} actual pair/modulus congruences")


def check_crt_conics() -> None:
    # All moduli are small enough for direct enumeration.
    for modulus in (3, 7, 9, 11, 21, 27, 33, 49, 63, 77):
        factors = factor(modulus)
        assert all(p % 4 == 3 for p in factors)
        expected = prod((p+1)*p**(a-1) for p, a in factors.items())
        assert conic_class_count(modulus, 1) == expected


def check_von_mangoldt_accounting() -> None:
    # Work prime-by-prime, so the identity is exact without floating logs.
    inert_primes = (3, 7, 11, 19, 23, 31)
    for m in range(2, 80):
        points, _ = build_points(m, 3)
        cofactors = [cofactor(points[i], points[j])
                     for i, j in combinations(range(m), 2)]
        for p in inert_primes:
            direct_exponent = 0
            depth_counts = 0
            for value in cofactors:
                exponent = 0
                x, y = value
                while x % p == 0 and y % p == 0:
                    exponent += 1
                    x //= p
                    y //= p
                direct_exponent += exponent
            power = p
            while power < m:
                count = sum(1 for i, j in combinations(range(m), 2)
                            if (i-j) % power == 0)
                assert count == e(m, power)
                depth_counts += count
                power *= p
            assert direct_exponent == depth_counts

    # Adding the 21-collision weight to its 3- and 7-collision weights
    # spends each of the two prime atoms twice on a pair divisible by 21.
    separate_atoms = {3: 1, 7: 1}
    naive_with_composite = {3: 2, 7: 2}
    assert naive_with_composite != separate_atoms


def check_joint_occupancy_bounds() -> None:
    moduli = (3, 7, 9, 11, 21, 27, 33, 49, 63, 77, 99)
    for m in range(2, 501):
        for modulus in moduli:
            factors = factor(modulus)
            q = prod((p+1)*p**(a-1) for p, a in factors.items())
            # The actual cyclic family has E(M,m) joint collisions and
            # respects the weaker universal Q(m)-class conic minimum.
            assert e(m, modulus) >= e(m, q)


def check_endpoint_scope() -> None:
    for m in range(4, 36):
        for n in range(1, 80):
            gap = 2*atan(2/(4*n*(n+1)+1))
            radius = prod(sqrt(4*(n+j)**2+1) for j in range(m))
            lower = (2/9) * 2**(m/2) * n**(m/2-2)
            assert gap >= 2/(9*n*n)
            assert gap*sqrt(radius) >= float(lower)


def main() -> None:
    check_actual_families()
    check_crt_conics()
    check_von_mangoldt_accounting()
    check_joint_occupancy_bounds()
    check_endpoint_scope()
    print("PASS: CRT conic counts and composite occupancy bounds.")
    print("PASS: exact prime-power accounting; composites add no von Mangoldt mass.")


if __name__ == "__main__":
    main()
