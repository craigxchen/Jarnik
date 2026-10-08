"""Exact checks of general endpoint contents and full Gaussian overlap costs."""

from itertools import combinations
from math import gcd, lcm, prod

from check_least_radius_formula import (
    conj, exact_div, gcd_all, gcd_gaussian, lcm_all, mul, norm,
)


def sub(a, b):
    return a[0] - b[0], a[1] - b[1]


def reduce_denominator(z):
    return exact_div(conj(z), gcd_gaussian(z, conj(z)))


def primitive_radius(rows):
    """Clear actual phase fractions, then independently remove tuple gcd."""
    nums, dens = [], []
    for z in rows:
        common = gcd_gaussian(z, conj(z))
        nums.append(exact_div(z, common))
        dens.append(exact_div(conj(z), common))
    den = lcm_all(dens)
    points = [mul(exact_div(den, b), a) for a, b in zip(nums, dens)]
    common = gcd_all(points)
    points = [exact_div(z, common) for z in points]
    assert all(norm(z) == norm(points[0]) for z in points)
    assert norm(gcd_all(points)) == 1
    return norm(points[0])


def all_edge_radius(rows):
    """Independent ordinary lcm using every relative phase, not an anchor star."""
    return lcm(*(norm(reduce_denominator(mul(a, conj(b))))
                 for a, b in combinations(rows, 2)))


def general_content():
    count = 0
    for x in range(1, 14):
        for y in range(1, 10):
            if gcd(x, y) != 1:
                continue
            for u in range(-12, 13):
                for v in range(-9, 10):
                    if gcd(u, v) != 1:
                        continue
                    d = y*u-x*v
                    c = gcd(d, v)
                    assert c == gcd(y, v) > 0
                    dp, vp = d//c, v//c
                    for D in range(1, 30):
                        h = gcd(dp, D)
                        expected = c*h*gcd(x*(dp//h)+(D//h)*vp, y)
                        actual = gcd(x*d+D*v, y*d)
                        assert actual == expected
                        count += 1
    return count


def cotangent_multiplier():
    count = 0
    for L in range(1, 18):
        for A in range(1, 26):
            candidates = [A+d for d in range(1, A*A+L*L+1)
                          if (A*A+L*L) % d == 0]
            for xs in combinations(candidates, 2):
                if (xs[0]*xs[1]+L*L) % (xs[1]-xs[0]):
                    continue
                g = gcd(A, L)
                x, y = A//g, L//g
                T = x*x+y*y
                ds = []
                for X in xs:
                    gi = gcd(X, L)
                    u, v = X//gi, L//gi
                    d = y*u-x*v
                    dp = d//gcd(d, v)
                    assert dp == (X-A)//gcd(g, gi)
                    ds.append(dp)
                E = lcm(*(d//gcd(d, T) for d in ds))
                assert g*g % E == 0
                assert T*E == lcm(T, *ds)
                for scale in (2, 5):
                    scaled = [(scale*X-scale*A)//gcd(scale*g, gcd(scale*X, scale*L))
                              for X in xs]
                    assert scaled == ds
                count += 1
    return count


def lcm_excess_fixture():
    # Q retains powers above the anchor prime; equation (5) still holds.
    pi = (2, 1)
    power = lambda n: prod_gaussian([pi]*n)
    B = power(2)
    As = [power(4), mul(power(3), (3, 2))]
    Qs = [exact_div(a, gcd_gaussian(B, a)) for a in As]
    assert norm(gcd_gaussian(B, Qs[0])) == 25
    assert norm(lcm_all([B]+As)) == norm(B)*norm(lcm_all(Qs))


def prod_gaussian(values):
    answer = (1, 0)
    for value in values:
        answer = mul(answer, value)
    return answer


def overlap_case(x, interiors, k):
    T, B = x*x+1, (x, 1)
    assert x > 0 and x % 2 == 0 and gcd(k, T) == 1
    Hs = [(d+x*v, v) for d, v in interiors]
    Cs = [(d+k*x*v, k*v) for d, v in interiors]
    Qs, Qps, ws = [], [], []
    for (d, v), H, C in zip(interiors, Hs, Cs):
        assert d > 0 and v > 0 and gcd(d, v) == 1
        # Conjugating every source denominator gives the B-oriented frame.
        A = conj(reduce_denominator(H))
        Ap = reduce_denominator(conj(C))
        assert norm(gcd_gaussian(B, A)) == gcd(d, T)
        assert norm(gcd_gaussian(B, Ap)) == gcd(d, T)
        Q = exact_div(A, gcd_gaussian(B, A))
        Qp = exact_div(Ap, gcd_gaussian(B, Ap))
        assert norm(Q) >= 2*x*v
        gamma = gcd(d, k)
        w = (k//gamma)*v
        assert gcd(d//gamma, w) == 1
        assert gcd(d//gamma, T) == gcd(d, T)
        assert norm(Qp) >= 2*x*w
        assert norm(Qp) <= (k*k if k % 2 else 2*k*k)*norm(Q)
        Qs.append(Q)
        Qps.append(Qp)
        ws.append(w)
    N = T*norm(lcm_all(Qs))
    Np = T*norm(lcm_all(Qps))
    assert N == primitive_radius([(1, 0), B]+Hs)
    D = k*T
    source = [(1, 0), B]+Hs
    target = [(x*u+(D-x*x)*v, u-x*v) for u, v in source]
    assert Np == primitive_radius(target)
    assert N == all_edge_radius(source)
    assert Np == all_edge_radius(target)
    omega_num, omega_den = prod(norm(q) for q in Qs), norm(lcm_all(Qs))
    assert omega_num % omega_den == 0
    omega = omega_num//omega_den
    # Algebraic part of (9), before the analytic Delta >= 1/x factor.
    assert omega*N >= T*(2*x)**len(interiors)
    target_num = prod(norm(q) for q in Qps)
    target_den = norm(lcm_all(Qps))
    assert target_num % target_den == 0
    target_omega = target_num//target_den
    # Exact algebraic input to (10); Delta >= 1/x is proved in the note.
    assert target_omega*Np == T*target_num
    assert target_omega*Np >= T*(2*x)**len(interiors)*prod(ws)
    factor = k*k if k % 2 else 2*k*k
    assert Np <= factor**len(interiors)*omega*N
    return omega


def overlap_checks():
    count, shared = 0, False
    for x in (2, 4, 6, 8, 12, 18):
        T = x*x+1
        options = [(d, v) for d in range(1, T+1) if T % d == 0
                   for v in range(1, 8) if gcd(d, v) == 1]
        for start in range(0, len(options), 3):
            rows = options[start:start+3]
            for k in (1, 2, 3, 4, 5):
                if gcd(k, T) != 1:
                    continue
                shared |= overlap_case(x, rows, k) > 1
                count += 1
    # T=325 includes 5^2; test source powers and an ordinary target content.
    overlap_case(18, [(1, 1), (1, 14), (25, 2), (325, 1)], 2)
    for k in (2, 3, 4, 7):
        overlap_case(18, [(1, 1), (125, 2), (650, 3), (k, 1)], k)
    assert shared
    return count+5


if __name__ == '__main__':
    print('general content cases:', general_content())
    print('integer-cotangent cliques:', cotangent_multiplier())
    lcm_excess_fixture()
    print('all-row overlap cases:', overlap_checks())
