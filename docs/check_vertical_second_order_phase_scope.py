"""Exact checks for vertical_second_order_phase_scope.md."""

from fractions import Fraction
from itertools import combinations, product
from math import gcd

from check_endpoint_vertical_line_gcd import exact_div, gcd_gaussian, norm


Gaussian = tuple[int, int]


def add(a: Gaussian, b: Gaussian) -> Gaussian:
    return a[0] + b[0], a[1] + b[1]


def mul(a: Gaussian, b: Gaussian) -> Gaussian:
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def ordinary_content(a: Gaussian) -> int:
    return gcd(abs(a[0]), abs(a[1]))


def gprod(values: list[Gaussian]) -> Gaussian:
    out = (1, 0)
    for value in values:
        out = mul(out, value)
    return out


def height_numerator(r: tuple[int, int, int, int]) -> int:
    pairs = sum(abs(r[i] + r[j]) for i, j in combinations(range(4), 2))
    s = sum(r)
    triples = sum(abs(s - r[j]) for j in range(4))
    return pairs + triples


def check_low_cost_lattice() -> int:
    count = 0
    for r in product(range(-2, 3), repeat=4):
        if not any(r) or height_numerator(r) >= 18:
            continue
        assert max(map(abs, r)) <= 2
        assert 125 * r[0] + 25 * r[1] + 5 * r[2] + r[3] != 0
        count += 1
    assert count == 380

    # Finite regression for the infinite triangle-inequality proof (14).
    for r in product(range(-5, 6), repeat=4):
        assert height_numerator(r) >= 6 * max(map(abs, r))
    return count


def low_cost_vectors() -> list[tuple[int, int, int, int]]:
    return [r for r in product(range(-2, 3), repeat=4)
            if any(r) and height_numerator(r) < 18]


def check_divided_differences() -> None:
    for x, p, q in ((11, 1, 1), (20, 4, 20), (37, 24, 100)):
        common = gcd(p, q)
        r = (q // common, -(p + q) // common, p // common, 0)
        ts = (x, x + p, x + p + q, x + p + q + 1)
        value = sum(Fraction(r[j], ts[j]) for j in range(3))
        assert value == Fraction(p * q * (p + q), common * ts[0] * ts[1] * ts[2])
        assert height_numerator(r) == 6 * (p + q) // common

    slopes = (1, 5, 25, 125)
    costs = []
    for i, j, k in combinations(range(4), 3):
        p, q = slopes[j] - slopes[i], slopes[k] - slopes[j]
        costs.append((p + q) // gcd(p, q))
    assert costs == [6, 31, 31, 6]

    pair_balance = (1, -1, -1, 1)
    assert height_numerator(pair_balance) == 8
    affine_nodes = (101, 107, 113, 119)
    assert sum(pair_balance) == 0
    assert sum(r * t for r, t in zip(pair_balance, affine_nodes)) == 0


def check_ap_content() -> None:
    for d, t, q in ((1, 20, 1), (3, 41, 7), (10, 100, 13)):
        ws = [(-d, t + j * q) for j in range(3)]
        u = mul(ws[0], ws[2])
        v = mul(ws[1], ws[1])
        raw = mul(u, (v[0], -v[1]))
        assert raw[1] == 2 * d * (t + q) * q * q
        common = gcd_gaussian(u, v)
        z = exact_div(raw, (norm(common), 0))
        assert z[1] != 0
        assert (2 * d * (t + q) * q * q) % norm(common) == 0
        content = ordinary_content(z)
        assert (2 * d * (t + q) * q * q) % (norm(common) * content) == 0


def check_endpoint_family() -> int:
    K = 3 * 7 * 13 * 31
    checked = 0
    for multiplier in (1, 3, 9, 21):
        M = K * multiplier
        assert M % 2 == 1
        ts = tuple(5**j * M for j in range(4))
        ws = [(-1, t) for t in ts]
        us = [exact_div(w, (1, 1)) for w in ws]
        norms = [norm(u) for u in us]
        assert all(value % 2 for value in norms)
        assert all(gcd(norms[i], norms[j]) == 1 for i, j in combinations(range(4), 2))
        assert all(norm(gcd_gaussian(us[i], us[j])) == 1
                   for i, j in combinations(range(4), 2))

        product_u = gprod(us)
        L = mul((1, 1), product_u)
        B = mul((1, 1), L)
        P = exact_div(B, (2, 0))
        assert P == mul((0, 1), product_u)

        vs = [exact_div(B, w) for w in ws]
        assert all(ordinary_content(v) == 1 for v in vs)
        qs = [add(P, v) for v in vs]
        assert all(norm(q) == norm(P) for q in qs)
        assert all(mul(vs[j], ws[j]) == B for j in range(4))

        common = P
        for qpoint in qs:
            common = gcd_gaussian(common, qpoint)
        assert norm(common) == 1

        # L is an exact lcm up to a unit: it is divisible by every w_j,
        # and its quotient has the expected norm product.
        for j, w in enumerate(ws):
            quotient = exact_div(L, w)
            expected = gprod([us[k] for k in range(4) if k != j])
            assert norm(quotient) == norm(expected)
        checked += 1
    return checked


def advance(z: Gaussian) -> Gaussian:
    x, y = z
    return 17 * x + 4 * y, 4 * x + y


def check_pell_endpoint() -> int:
    # Coefficient pairs (rational part, coefficient of s) for
    # CE, BE, AE, BC, where s^2+4s-1=0.
    leading = ((5, 10), (5, -5), (1, -2), (-5, 45))
    vectors = low_cost_vectors()
    for r in vectors:
        rational = sum(r[j] * leading[j][0] for j in range(4))
        irrational = sum(r[j] * leading[j][1] for j in range(4))
        assert (rational, irrational) != (0, 0)

    a, b, c, e = (0, 1), (1, 3), (4, -3), (2, 1)
    g = (1, 2)
    checked = 0
    for n in range(26):
        if n and n % 5 == 0:
            ws = [mul((2, 0), mul(g, mul(x, y)))
                  for x, y in ((a, b), (a, c), (b, c), (a, e))]
            ts = tuple(w[1] for w in ws)
            assert all(t > 0 for t in ts)
            assert len(set(ts)) == 4
            assert all(sum(Fraction(r[j], ts[j]) for j in range(4)) != 0
                       for r in vectors)
            checked += 1
        a, b, c, e = map(advance, (a, b, c, e))
    return checked


if __name__ == "__main__":
    low_cost = check_low_cost_lattice()
    check_divided_differences()
    check_ap_content()
    families = check_endpoint_family()
    pell = check_pell_endpoint()
    print(f"PASS: {low_cost} low-cost vectors; divided differences and AP content")
    print(f"PASS: {families} primitive endpoint instances with exact Gaussian identities")
    print(f"PASS: {pell} bounded-arc Pell instances and all limiting low-cost phases")
