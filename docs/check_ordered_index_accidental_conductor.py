"""Exact checks for accidental conductor in the ordered five-row index."""

from functools import reduce
from itertools import combinations, product
from math import gcd

from check_joint_affine_relation_lattice import (
    ggcd, gnorm, gsub, primitive_points,
)
from check_ordered_affine_circuit_conductor_index import (
    band_index, gvaluation, valuation,
)


def vp(n, p, cap=10**9):
    if n == 0:
        return cap
    answer = 0
    while n % p == 0:
        n //= p
        answer += 1
    return answer


def triangle_value(bits, depths, triple, e):
    values = [bits[i] for i in triple]
    base = e if len(set(values)) == 1 else 0
    return base + sum(depths[i, j] for i, j in combinations(triple, 2)
                      if bits[i] == bits[j])


def gcd_value(bits, depths, rows, e):
    return min(triangle_value(bits, depths, triple, e)
               for triple in combinations(rows, 3))


def mu(bits, depths, rows):
    return min(depths[i, j] for i, j in combinations(rows, 2)
               if bits[i] == bits[j])


def nu(depths, rows):
    return min(sum(depths[i, j] for i, j in combinations(triple, 2))
               for triple in combinations(rows, 3))


def gap(bits, e):
    endpoints = [e*bits[i] for i in (0, 4)]
    interior = [e*bits[i] for i in (1, 2, 3)]
    return max(0, min(endpoints)-max(interior),
               min(interior)-max(endpoints))


def unit_depths(units, bits, p):
    answer = {}
    for i, j in combinations(range(5), 2):
        if bits[i] == bits[j]:
            answer[i, j] = vp(units[i]-units[j], p)
    return answer


def exhaustive_two_level():
    checked = accidental = 0
    # Actual integer units modulo powers of five automatically give
    # ultrametric equal-level depth tables.
    unit_choices = (1, 2, 3, 6, 11)
    for bits in product((0, 1), repeat=5):
        if len(set(bits)) == 1:
            continue
        for units in product(unit_choices, repeat=5):
            if any(units[i] == units[j] for i, j in combinations(range(5), 2)
                   if bits[i] == bits[j]):
                continue
            depths = unit_depths(units, bits, 5)
            e = 3
            full = gcd_value(bits, depths, range(5), e)
            left = gcd_value(bits, depths, range(4), e)
            right = gcd_value(bits, depths, range(1, 5), e)
            middle = triangle_value(bits, depths, (1, 2, 3), e)
            q = full + middle - left - right
            assert q >= gap(bits, e)

            interior = bits[1:4]
            if len(set(interior)) == 2:
                repeated = next(pair for pair in combinations((1, 2, 3), 2)
                                if bits[pair[0]] == bits[pair[1]])
                predicted = (mu(bits, depths, range(5))+depths[repeated]
                             -mu(bits, depths, range(4))
                             -mu(bits, depths, range(1, 5)))
                assert q == predicted
                assert gap(bits, e) == 0
                accidental += q > 0
            elif gap(bits, e) == 0:
                # The only primitive case is a singleton endpoint.
                s = sum(depths[i, j] for i, j in combinations((1, 2, 3), 2))
                m123 = min(depths[i, j]
                           for i, j in combinations((1, 2, 3), 2))
                if bits[0] != bits[1]:
                    predicted = (mu(bits, depths, range(1, 5))+s-m123
                                 -nu(depths, range(1, 5)))
                else:
                    assert bits[4] != bits[1]
                    predicted = (mu(bits, depths, range(4))+s
                                 -nu(depths, range(4))-m123)
                assert q == predicted
                accidental += q > 0
            checked += 1
    assert accidental
    return checked, accidental


def local_unbounded():
    p, e = 5, 4
    bits = (0, 1, 1, 0, 0)
    for k in range(1, 9):
        precision = e+k+6
        modulus = p**precision
        units = [1, 1, 1+p**k, 2, 3]
        xs, ys = [], []
        for bit, unit in zip(bits, units):
            t = e*bit
            inverse = pow(unit, -1, modulus)
            x = (p**t*unit) % modulus
            y = (p**(e-t)*inverse) % modulus
            assert (x*y-p**e) % modulus == 0
            xs.append(x)
            ys.append(y)
        depths = unit_depths(units, bits, p)
        assert depths[1, 2] == k
        assert depths[0, 3] == depths[0, 4] == depths[3, 4] == 0
        # Check both split-coordinate difference valuations.
        for i, j in combinations(range(5), 2):
            ti, tj = e*bits[i], e*bits[j]
            expected_x = min(ti, tj) if ti != tj else ti+depths[i, j]
            expected_y = min(e-ti, e-tj) if ti != tj else e-ti+depths[i, j]
            assert vp((xs[i]-xs[j]) % modulus, p, precision) == expected_x
            assert vp((ys[i]-ys[j]) % modulus, p, precision) == expected_y
        for triple in combinations(range(5), 3):
            i, j, ell = triple
            determinant = ((xs[j]-xs[i])*(ys[ell]-ys[i])
                           -(xs[ell]-xs[i])*(ys[j]-ys[i])) % modulus
            assert vp(determinant, p, precision) == triangle_value(
                bits, depths, triple, e)
        full = gcd_value(bits, depths, range(5), e)
        left = gcd_value(bits, depths, range(4), e)
        right = gcd_value(bits, depths, range(1, 5), e)
        middle = triangle_value(bits, depths, (1, 2, 3), e)
        assert gap(bits, e) == 0
        assert full+middle-left-right == k


def forced_factor(points):
    n = gnorm(points[0])
    ao = gnorm(ggcd(points[0], points[4]))
    ai = gnorm(reduce(ggcd, points[1:4]))
    x = ao*ai
    return x//gcd(x, n), ao, ai


def global_fixtures():
    points = [(446,23),(439,82),(343,286),(329,302),(302,329)]
    n = 199445
    assert all(gnorm(z) == n for z in points)
    assert gnorm(reduce(ggcd, points)) == 1
    q, g, local = band_index(points, all_minors=True)
    assert (q, g, local, gcd(q, n)) == (220, 6, [6, 6], 5)
    pi = (1, 2)
    t = [gvaluation(z, pi) for z in points]
    assert t == [0, 0, 1, 0, 1]
    depths = {(i, j):gvaluation(gsub(points[i], points[j]), pi)-t[i]
              for i, j in combinations(range(5), 2) if t[i] == t[j]}
    assert depths == {(0,1):0, (0,3):0, (1,3):1, (2,4):0}
    assert valuation(q, 5) == 1
    assert forced_factor(points) == (1, 353, 1)
    chord2 = gnorm(gsub(points[4], points[0]))
    assert chord2 == 114372 and 16*chord2**2 > n

    quartic = [primitive_points(12)[j] for j in (3, 0, 4, 1, 2)]
    n = gnorm(quartic[0])
    q, _, _ = band_index(quartic, all_minors=True)
    assert q == 2937 and gcd(q, n) == 89
    assert valuation(n, 89) == 2 and valuation(q, 89) == 1
    assert [gvaluation(z, (5,8)) for z in quartic] == [0,1,1,1,2]
    assert forced_factor(quartic)[0] == 1
    chord2 = gnorm(gsub(quartic[4],quartic[0]))
    assert 16*chord2**2 > n


def main():
    checked, accidental = exhaustive_two_level()
    local_unbounded()
    global_fixtures()
    print(f"PASS: {checked} two-level tables, {accidental} accidental cases, "
          "local unbounded family, and two global fixtures")


if __name__ == "__main__":
    main()
