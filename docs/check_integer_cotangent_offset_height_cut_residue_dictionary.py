"""Literal cut/residue formulas for normalized cotangent offset height.

The Gaussian source tuple is constructed independently of the integer
offset calculation.  Tests retain split orientations, inert content,
ramified parity, common scaling, and the formal full-fair count.
"""

from itertools import combinations, product
from math import atan, ceil, gcd, isqrt, lcm, log

from check_integer_cotangent_normalization import primitive_tuple
from check_least_radius_formula import conj, divmod_gaussian, exact_div, mul, norm


def factor(n: int) -> dict[int, int]:
    n = abs(n)
    result: dict[int, int] = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            result[p] = result.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        result[n] = result.get(n, 0) + 1
    return result


def vp(n: int, p: int) -> int:
    assert n != 0
    return factor(n).get(p, 0)


def edge(q: int, length: int) -> tuple[int, int, int, int]:
    g = gcd(abs(q), length)
    r = length // g
    a = q // g
    eta = int(bool(a % 2 and r % 2))
    n = (a * a + r * r) // (2**eta)
    assert n % 2
    numerator = (a, r)
    h = exact_div(numerator, (1, 1)) if eta else numerator
    assert norm(h) == n
    assert r == (h[0] + h[1] if eta else h[1])
    unit = (0, 1) if eta else (1, 0)
    assert mul(numerator, conj(h)) == mul(mul(unit, h), conj(numerator))
    return g, r, eta, n


def gaussian_prime(p: int) -> tuple[int, int]:
    assert p % 4 == 1
    for a in range(1, isqrt(p) + 1):
        b2 = p - a * a
        b = isqrt(b2)
        if b > 0 and b * b == b2:
            return a, b
    raise AssertionError(p)


def gaussian_v(z: tuple[int, int], pi: tuple[int, int]) -> int:
    value = z
    count = 0
    while True:
        quotient, remainder = divmod_gaussian(value, pi)
        if remainder != (0, 0):
            return count
        count += 1
        value = quotient


def check_clique(xs: tuple[int, ...], length: int) -> tuple[int, int]:
    assert length > 0 and len(xs) >= 3
    assert all(a < b for a, b in zip(xs, xs[1:]))
    a = xs[0]
    q_anchor = []
    for x in xs[1:]:
        numerator = a * x + length * length
        denominator = a - x
        assert numerator % denominator == 0
        q_anchor.append(numerator // denominator)
    for x, y in combinations(xs[1:], 2):
        assert (x * y + length * length) % (x - y) == 0

    ds = tuple(x - a for x in xs[1:])
    h = lcm(*ds) // gcd(*ds)
    s = a * a + length * length
    assert all(s % d == 0 for d in ds)
    cs = tuple(s // d for d in ds)
    assert h == s // (gcd(*ds) * gcd(*cs))

    q_all = [*xs]
    q_all.extend((x * y + length * length) // (x - y)
                 for x, y in combinations(xs, 2))
    n = lcm(*(edge(q, length)[3] for q in q_all))
    rows = primitive_tuple(xs, length)
    assert all(norm(z) == n for z in rows)
    ps = set().union(*(factor(v) for v in
                       (h, n, length, *ds, *xs, *q_all)))
    assert all(p % 4 == 1 for p in factor(n))

    direct_core = 0.0
    gamma = 0.0
    for p in sorted(ps):
        if p in factor(n):
            pi = gaussian_prime(p)
            es = [gaussian_v(z, pi) for z in rows]
            assert min(es) == 0 and max(es) == vp(n, p)
        else:
            es = [0] * len(rows)
        tvals = []
        corrected = []
        for pos, (x, q, d) in enumerate(zip(xs[1:], q_anchor, ds), start=2):
            ex, rx, etax, nx = edge(x, length)
            ea, ra, etaa, na = edge(a, length)
            eq, rq, etaq, nq = edge(q, length)
            if p in factor(n):
                t = (abs(es[pos] - es[0])
                     + abs(es[1] - es[0])
                     - abs(es[pos] - es[1])) // 2
                assert 2 * t == (vp(nx, p) + vp(na, p) - vp(nq, p))
            else:
                t = 0
                assert vp(nx, p) == vp(na, p) == vp(nq, p) == 0
            delta = (etax + etaa - etaq) // 2 if p == 2 else 0
            assert delta in (0, 1)
            if p == 2:
                assert etax + etaa - etaq == 2 * delta
            expression = t + vp(ea, p) + vp(rq, p) - vp(rx, p) + delta
            assert expression == vp(d, p), (xs, length, p, pos, expression)
            tvals.append(t)
            corrected.append(t + vp(rq, p) - vp(rx, p) + delta)
        assert max(corrected) - min(corrected) == vp(h, p)
        if p in factor(n):
            core_p = max(tvals) - min(tvals)
            direct_core += core_p * log(p)
            ei = es[2:]
            lost = abs(es[1] - es[0]) - core_p
            removed = max(es) - min(es) - (max(ei) - min(ei))
            assert 0 <= lost <= removed
            gamma += removed * log(p)

    assert abs(log(h) - direct_core) <= 2 * log(length) + log(2) + 1e-10
    assert 0 <= log(edge(a, length)[3]) - direct_core <= gamma + 1e-10
    cstar = 2 * n**0.25 * atan(length / a)
    ell = ceil((len(xs) - 1) / (4 * ceil(cstar / 2)))
    assert gamma <= log(n) / ell + 1e-10
    lower = (0.5 - 1 / ell) * log(n) - 2 * log(length) - 2 * log(cstar)
    assert log(h) + 1e-10 >= lower
    return h, n


def check_formal_fair() -> int:
    count = 0
    for total_rows in range(4, 10):
        cuts = [(0, *bits) for bits in product((0, 1), repeat=total_rows - 1)
                if any(bits)]
        assert len(cuts) == 2 ** (total_rows - 1) - 1
        good = sum(cut[1] == 1 and len(set(cut[2:])) == 2
                   for cut in cuts)
        assert good == 2 ** (total_rows - 2) - 2
        count += len(cuts)
    return count


def main() -> None:
    fixtures = (
        ((85, 182, 210), 70, 12125, 485),
        ((157, 182, 447), 1, 290, 212298125),
        ((3, 4, 9), 3, 6, 25),
        ((4, 9, 12), 12, 40, 25),
    )
    checks = 0
    for xs, length, h, n in fixtures:
        assert check_clique(xs, length) == (h, n)
        checks += 1
        for scale in (2, 3, 7):
            assert check_clique(tuple(scale * x for x in xs),
                                scale * length) == (h, n)
            checks += 1
    assert 3 in factor(fixtures[2][2]) and 3 not in factor(fixtures[2][3])
    assert 2 in factor(fixtures[2][2]) and 2 not in factor(fixtures[2][3])
    searched = 0
    for length in range(1, 9):
        for a in range(length, 31):
            s = a * a + length * length
            ds = [d for d in range(1, 61 - a) if s % d == 0]
            for d1, d2 in combinations(ds, 2):
                x, y = a + d1, a + d2
                if (x * y + length * length) % (x - y) == 0:
                    check_clique((a, x, y), length)
                    searched += 1
    assert searched == 682
    cuts = check_formal_fair()
    print(f"PASS: {checks} fixtures, {searched} searched cliques, {cuts} fair cuts")


if __name__ == "__main__":
    main()
