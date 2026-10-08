"""Exact checks for the oriented opposite-twist descent of the gap K."""

from functools import reduce
from itertools import product
from math import gcd
from random import Random

from check_ordered_affine_circuit_conductor_index import (
    gconj,
    gdiv_exact,
    ggcd,
    gmul,
    gnorm,
    gpower,
)


O = (0, 4)
I = (1, 2, 3)


def interval(t, rows):
    vals = [t[j] for j in rows]
    return min(vals), max(vals)


def oriented_gap(t):
    """Return (h, side): side=1 means outer above, -1 interior above."""
    a, b = interval(t, O)
    u, v = interval(t, I)
    if a > v:
        return a - v, 1
    if u > b:
        return u - b, -1
    return 0, 0


def descended_allocations(t, e):
    h, side = oriented_gap(t)
    if side == 1:
        # alpha contains pi^h: endpoint pi-valuations decrease.
        new = [t[j] - h if j in O else t[j] for j in range(5)]
    elif side == -1:
        # alpha contains bar(pi)^h, so bar(alpha) contains pi^h:
        # interior pi-valuations decrease.
        new = [t[j] - h if j in I else t[j] for j in range(5)]
    else:
        new = list(t)
    return h, side, new, e - h


def check_allocation(t, e):
    h, _, new, enew = descended_allocations(t, e)
    assert min(new) == 0 and max(new) == enew
    assert all(0 <= q <= enew for q in new)
    assert oriented_gap(new)[0] == 0
    return h


def strict_minor_order(points):
    """Exact sufficient test: every later ray is CCW by less than pi/2."""
    for i, z in enumerate(points):
        for w in points[i + 1 :]:
            cross = z[0] * w[1] - z[1] * w[0]
            dot = z[0] * w[0] + z[1] * w[1]
            assert cross > 0 and dot > 0


def check_constructed(points, alpha, k):
    assert gnorm(alpha) == k
    a, b = alpha
    big_a, big_b = a * a - b * b, 2 * a * b
    assert big_a * big_a + big_b * big_b == k * k
    assert gcd(abs(big_a), abs(big_b)) == 1
    descended = []
    for j, z in enumerate(points):
        divisor = alpha if j in O else gconj(alpha)
        descended.append(gdiv_exact(z, divisor))
    n = gnorm(points[0])
    assert all(gnorm(z) == n for z in points)
    assert all(gnorm(z) == n // k for z in descended)
    assert reduce(ggcd, descended) in ((1, 0), (-1, 0), (0, 1), (0, -1))
    return descended


def main():
    allocations = 0
    nonzero_gaps = 0
    for e in range(1, 9):
        for t in product(range(e + 1), repeat=5):
            if min(t) != 0 or max(t) != e:
                continue
            nonzero_gaps += check_allocation(t, e) > 0
            allocations += 1

    rng = Random(922603)
    pis = ((2, 1), (3, 2), (4, 1))  # norms 5, 13, 17
    units = ((1, 0), (0, 1), (-1, 0), (0, -1))
    actual = 0
    actual_nontrivial = 0
    for _ in range(800):
        points = [rng.choice(units) for _ in range(5)]
        alpha = (1, 0)
        k = 1
        for pi in pis:
            e = rng.randrange(1, 6)
            t = [rng.randrange(e + 1) for _ in range(5)]
            # Remove the full common layer first; this is the primitive
            # allocation seen by the theorem at this prime.
            lo, hi = min(t), max(t)
            t = [q - lo for q in t]
            e = hi - lo
            if e == 0:
                continue
            h, side, _, _ = descended_allocations(t, e)
            factor = gpower(pi if side == 1 else gconj(pi), h)
            alpha = gmul(alpha, factor)
            k *= gnorm(pi) ** h
            for j in range(5):
                points[j] = gmul(
                    points[j],
                    gmul(gpower(pi, t[j]), gpower(gconj(pi), e - t[j])),
                )
        check_constructed(points, alpha, k)
        actual_nontrivial += k > 1
        actual += 1

    fixture = [(4, -33), (9, -32), (12, -31), (23, -24), (24, -23)]
    strict_minor_order(fixture)
    descended = check_constructed(fixture, (-1, 2), 5)
    assert descended == [(-14, 5), (11, 10), (10, 11), (5, 14), (-14, -5)]

    # Boundary numerators in (9) are literal nonzero Gaussian integers.
    alpha = (-1, 2)
    a, b = alpha
    big_a, big_b = a * a - b * b, 2 * a * b
    for x in (descended[0], descended[4]):
        for y in descended[1:4]:
            c, d = gmul(y, gconj(x))
            original_x = gmul(alpha, x)
            original_y = gmul(gconj(alpha), y)
            cross = original_x[0] * original_y[1] - original_x[1] * original_y[0]
            assert cross == big_a * d - big_b * c
            numerator = (
                alpha[0] * x[0] - alpha[1] * x[1]
                - (gconj(alpha)[0] * y[0] - gconj(alpha)[1] * y[1]),
                alpha[0] * x[1] + alpha[1] * x[0]
                - (gconj(alpha)[0] * y[1] + gconj(alpha)[1] * y[0]),
            )
            assert numerator != (0, 0)

    assert nonzero_gaps and actual_nontrivial
    print(
        f"PASS: {allocations} primitive allocations ({nonzero_gaps} gaps), "
        f"{actual} random Gaussian tuples, and the ordered K=5 fixture"
    )


if __name__ == "__main__":
    main()
