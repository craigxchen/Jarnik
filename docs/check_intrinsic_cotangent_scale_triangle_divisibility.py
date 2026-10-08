"""Exact intrinsic cotangent-scale divisibilities and inert-prime sharpness."""

from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import comb, gcd, isqrt, lcm, prod

from check_least_radius_formula import (
    conj, exact_div, gcd_all, gcd_gaussian, lcm_all, mul, norm,
)


def sub(a, b):
    return a[0]-b[0], a[1]-b[1]


def dot(a, b):
    return a[0]*b[0]+a[1]*b[1]


def det(a, b):
    return a[0]*b[1]-a[1]*b[0]


def valuation(value, p):
    assert value
    count = 0
    while value % p == 0:
        value //= p
        count += 1
    return count


def primitive_phases(half_angles):
    nums, dens = [], []
    for H in half_angles:
        common = gcd_gaussian(H, conj(H))
        nums.append(exact_div(H, common))
        dens.append(exact_div(conj(H), common))
    denominator = lcm_all(dens)
    points = [mul(exact_div(denominator, d), n) for n, d in zip(nums, dens)]
    common = gcd_all(points)
    return [exact_div(z, common) for z in points]


def all_edge_radius(half_angles):
    edges = []
    for a, b in combinations(half_angles, 2):
        H = mul(a, conj(b))
        common = gcd_gaussian(H, conj(H))
        edges.append(norm(exact_div(conj(H), common)))
    return lcm(*edges)


def check_tuple(points):
    assert len(points) >= 3 and len(set(points)) == len(points)
    N = norm(points[0])
    assert N % 2 == 1
    assert all(norm(z) == N for z in points)
    assert norm(gcd_all(points)) == 1
    pairs = list(combinations(range(len(points)), 2))
    triples = list(combinations(range(len(points)), 3))
    cot, den, halves, tri_lcms, tri_contents = {}, {}, {}, {}, {}
    for i, j in pairs:
        cross = det(points[i], points[j])
        if cross == 0:
            assert points[j] == (-points[i][0], -points[i][1])
            cot[i, j] = Fraction(0)
        else:
            cot[i, j] = Fraction(N+dot(points[i], points[j]), cross)
        den[i, j] = cot[i, j].denominator
        a, s = cot[i, j].numerator, den[i, j]
        eps = 2 if a % 2 and s % 2 else 1
        pair_norm = N//norm(gcd_gaussian(points[i], points[j]))
        assert (a*a+s*s)//eps == pair_norm
        assert eps == (2 if points[i][0] % 2 != points[j][0] % 2 else 1)
        assert norm(sub(points[i], points[j]))*eps*pair_norm == 4*N*s*s
        cofactor = exact_div(sub(points[i], points[j]), gcd_gaussian(points[i], points[j]))
        assert norm(cofactor)*eps == 4*s*s
    for i, j, k in triples:
        D = det(sub(points[j], points[i]), sub(points[k], points[i]))
        assert D and D % 2 == 0
        halves[i, j, k] = abs(D)//2
        tri_lcms[i, j, k] = lcm(den[i, j], den[i, k], den[j, k])
        assert halves[i, j, k] % tri_lcms[i, j, k] == 0
        tri_contents[i, j, k] = norm(gcd_all([points[i], points[j], points[k]]))
        eta = 2 if len({points[t][0] % 2 for t in (i, j, k)}) == 1 else 1
        assert halves[i, j, k] == eta*tri_contents[i, j, k]*den[i, j]*den[i, k]*den[j, k]
    for i, j in pairs:
        incident = []
        for k in range(len(points)):
            if k in (i, j):
                continue
            a, b = sub(points[j], points[i]), sub(points[k], points[i])
            D, B = det(a, b), norm(b)-dot(a, b)
            assert D*cot[i, j] == B
            assert B % 2 == 0
            assert (D//2) % den[i, j] == 0
            incident.append(abs(D)//2)
        assert reduce(gcd, incident) % den[i, j] == 0
    L = lcm(*den.values())
    middle, right = prod(tri_lcms.values()), prod(halves.values())
    assert middle % (L**(len(points)-2)) == 0
    assert right % middle == 0
    m = len(points)
    bm = comb(m//2, 3)+comb((m+1)//2, 3)
    assert prod(tri_contents.values()) % (N**bm) == 0
    assert right % ((2*N)**bm*prod(den.values())**(m-2)) == 0
    return L, den, halves, tri_lcms


def small_actual_circles():
    count = 0
    for N in range(1, 180, 2):
        bound = isqrt(N)
        rows = [(x, y) for x in range(-bound, bound+1)
                for y in range(-bound, bound+1) if x*x+y*y == N]
        if not rows or norm(gcd_all(rows)) != 1:
            continue
        # Full circles deliberately test antipodes as well as short subsets.
        check_tuple(rows)
        count += 1
        positive = [z for z in rows if z[0] > 0 and z[1] >= 0]
        if len(positive) >= 3:
            for subset in combinations(positive, min(5, len(positive))):
                common = gcd_all(subset)
                check_tuple([exact_div(z, common) for z in subset])
                count += 1
    return count


def sharpness():
    count = 0
    for m in range(3, 10):
        p = 7 if m < 7 else 11
        assert p > m and p % 4 == 3
        for h in range(1, 5):
            ts = [0, p**h]+list(range(1, m-1))
            for H in (1, p**h+1, p**(h+2)+1):
                assert gcd(H, p) == 1
                half_angles = [(H, t) for t in ts]
                points = primitive_phases(half_angles)
                N = norm(points[0])
                assert N == all_edge_radius(half_angles)
                assert N % p
                L, den, halves, tri_lcms = check_tuple(points)
                assert valuation(L, p) == h
                for (i, j), s in den.items():
                    assert Fraction(H*H+ts[i]*ts[j], H*(ts[j]-ts[i])).denominator == s
                    assert valuation(s, p) == (h if (i, j) == (0, 1) else 0)
                for triple, a in halves.items():
                    expected = h if 0 in triple and 1 in triple else 0
                    assert valuation(a, p) == expected
                    assert valuation(tri_lcms[triple], p) == expected
                assert valuation(prod(halves.values()), p) == (m-2)*h
                count += 1
    return count


if __name__ == '__main__':
    for m in range(3, 101):
        bm = comb(m//2, 3)+comb((m+1)//2, 3)
        assert comb(m, 3)-4*bm == (m-2)*(m//2)
    print('actual odd-norm circle fixtures:', small_actual_circles())
    print('fully normalized inert-prime sharpness fixtures:', sharpness())
