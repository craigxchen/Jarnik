"""Exact checks for vertical_lcm_pell_counterfamily.md."""

from fractions import Fraction
from math import gcd

from check_endpoint_vertical_line_gcd import (
    add, exact_div, gcd_gaussian, gcd_many, lcm_many, mul, norm, sub,
)


def advance(z):
    x, y = z
    return 17 * x + 4 * y, 4 * x + y


def content(z):
    return gcd(abs(z[0]), abs(z[1]))


def check(indices=40):
    a, b, c = (0, 1), (1, 3), (4, -3)
    e = (2, 1)
    g, h = (1, 2), (5, 5)
    histories = []
    ratios = []
    for n in range(indices):
        assert all(content(v) == 1 for v in (a, b, c))
        assert a[0] * b[1] - a[1] * b[0] == -1
        assert norm(gcd_gaussian(a, b)) == 1
        p = mul(g, mul(a, mul(b, c)))
        vs = [c, b, a]
        ws = [mul((2, 0), mul(g, mul(x, y)))
              for x, y in ((a, b), (a, c), (b, c))]
        points = [p] + [add(p, mul((10, 0), v)) for v in vs]
        assert all(norm(q) == norm(p) for q in points)
        assert len(set(points)) == 4
        assert [content(sub(q, p)) for q in points[1:]] == [10, 10, 10]
        assert all(w[0] == -10 for w in ws)
        assert all(mul(v, w) == mul((2, 0), p) for v, w in zip(vs, ws))
        assert norm(lcm_many(ws)) == 4 * norm(p)
        assert norm(gcd_many(points)) == norm(h)
        for q in points:
            exact_div(q, h)

        reduced = [exact_div(q, h) for q in points]
        assert norm(gcd_many(reduced)) == 1
        p_reduced = reduced[0]
        chords = [sub(q, p_reduced) for q in reduced[1:]]
        ds = [content(v) for v in chords]
        assert ds == [1, 2, 1]
        primitive = [(v[0] // d, v[1] // d) for v, d in zip(chords, ds)]
        assert norm(gcd_many(primitive)) == 1
        reduced_ws = [exact_div(mul((2, 0), p_reduced), v) for v in primitive]
        assert norm(lcm_many(reduced_ws)) == 4 * norm(p_reduced)
        assert all(w[0] == -d for w, d in zip(reduced_ws, ds))
        assert [Fraction(w[1], d) for w, d in zip(reduced_ws, ds)] == [
            Fraction(w[1], 10) for w in ws
        ]

        ts = [w[1] for w in ws]
        histories.append(ts)
        if n >= 2:
            assert ts == [322 * x - y for x, y in zip(histories[-2], histories[-3])]
        if n >= 1:
            assert 0 < ts[0] < ts[1] < ts[2]
            assert all(t * t <= 4 * norm(p) for t in ts)
            # Square of L/(min(t)/d)^2, computed without floating point.
            ratios.append(Fraction(4 * norm(p) * 10**4, min(ts)**4))
            if len(ratios) >= 2:
                assert ratios[-1] < ratios[-2]
        if n % 5 == 0:
            four_complements = [mul(c, e), mul(b, e), mul(a, e), mul(b, c)]
            assert all(content(v) == 1 for v in four_complements)
            four_ws = ws + [mul((2, 0), mul(g, mul(a, e)))]
            four_a = mul(mul((2, 0), p), e)
            four_p = mul(p, e)
            assert all(w[0] == -10 for w in four_ws)
            assert all(mul(v, w) == four_a
                       for v, w in zip(four_complements, four_ws))
            five_points = [four_p] + [add(four_p, mul((10, 0), v))
                                      for v in four_complements]
            assert len(set(five_points)) == 5
            assert all(norm(q) == norm(four_p) for q in five_points)
            assert norm(gcd_many(four_complements)) <= 50
            four_lcm_norm = norm(lcm_many(four_ws))
            assert norm(four_a) <= 50 * four_lcm_norm <= 50 * norm(four_a)
            if n:
                assert len({w[1] for w in four_ws}) == 4
                assert all(w[1] > 0 for w in four_ws)
        a, b, c = map(advance, (a, b, c))
        e = advance(e)
    print(f"Verified {indices} exact indices: equal content 10; primitive contents (1,2,1).")
    print("Exact lcm reconstruction and Gaussian gcd 5(1+i) verified at every index.")
    print(f"Squared quadratic-bound ratio: {float(ratios[0]):.6g} -> {float(ratios[-1]):.6g}.")
    print("Four-cofactor quadratic-growth variant verified on indices divisible by five.")


if __name__ == "__main__":
    check()
