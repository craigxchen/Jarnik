"""Exact checks for the four-finite-cotangent Vandermonde alternative."""

from itertools import combinations
from math import gcd

from check_integer_cotangent_normalization import (
    check_clique,
    edge_norm,
    edge_quotient,
)


def is_clique(xs, scale):
    return all(
        (x * y + scale * scale) % (x - y) == 0
        for x, y in combinations(xs, 2)
    )


def audit(xs, scale):
    assert len(xs) == 4 and tuple(sorted(xs)) == xs and xs[0] > 0
    A = xs[0]
    S = A * A + scale * scale
    ds = [x - A for x in xs[1:]]
    assert all(S % d == 0 for d in ds)
    cs = [S // d for d in ds]
    es = [2 * A + d + c for d, c in zip(ds, cs)]
    assert all(d * e == x * x + scale * scale for d, e, x in zip(ds, es, xs[1:]))

    rs = {}
    gs = {}
    for i, j in combinations(range(3), 2):
        g = gcd(ds[i], ds[j])
        r = abs(ds[i] - ds[j]) // g
        gs[i, j] = g
        rs[i, j] = r
        assert es[i] % r == 0 and es[j] % r == 0
        assert (es[i] - es[j]) % r == 0

        hi, lo = (i, j) if ds[i] > ds[j] else (j, i)
        q = edge_quotient(xs[hi + 1], xs[lo + 1], scale)
        assert q > xs[lo + 1]
        assert (ds[hi] // g) * (es[hi] // r) == q + xs[hi + 1]
        assert (ds[lo] // g) * (es[lo] // r) == q - xs[lo + 1]

    rprod = 1
    gprod = 1
    vander = 1
    for i, j in combinations(range(3), 2):
        rprod *= rs[i, j]
        gprod *= gs[i, j]
        vander *= abs(ds[i] - ds[j])
    eprod = es[0] * es[1] * es[2]
    dprod = ds[0] * ds[1] * ds[2]
    assert eprod % rprod == 0
    invariant = eprod // rprod
    assert invariant > 0
    assert invariant * rprod == eprod

    # Complement symmetry preserves normalized gaps and the invariant.
    dual_ds = cs
    dual_es = [2 * A + c + d for c, d in zip(cs, ds)]
    dual_rprod = 1
    for i, j in combinations(range(3), 2):
        dual_g = gcd(dual_ds[i], dual_ds[j])
        dual_r = abs(dual_ds[i] - dual_ds[j]) // dual_g
        assert dual_r == rs[i, j]
        dual_rprod *= dual_r
    assert dual_es == es and eprod // dual_rprod == invariant

    n = check_clique(xs, scale)
    edge_values = list(xs)
    edge_values.extend(
        edge_quotient(x, y, scale) for x, y in combinations(xs, 2)
    )
    assert all(n % edge_norm(q, scale) == 0 for q in edge_values)

    lhs = 1
    for x in xs[1:]:
        lhs *= x * x + scale * scale
    assert lhs == dprod * eprod
    # Squared exact form of lhs <= 8 L^3 N^(3/2) V.
    assert lhs * lhs <= 64 * scale**6 * n**3 * vander**2
    return invariant


def main():
    checked = 0
    invariants = set()
    for scale in range(1, 13):
        values = range(1, 61)
        for xs in combinations(values, 4):
            if not is_clique(xs, scale):
                continue
            invariants.add(audit(xs, scale))
            checked += 1

    progression = 0
    for t in range(1, 101):
        scale = 6
        xs = tuple(scale * (t + i) for i in range(4))
        assert is_clique(xs, scale)
        audit(xs, scale)
        d = [x - xs[0] for x in xs[1:]]
        vander = (d[1] - d[0]) * (d[2] - d[0]) * (d[2] - d[1])
        assert vander == 2 * scale**3
        progression += 1

    print(f"PASS: {checked} positive four-cotangent cliques.")
    print(f"PASS: {len(invariants)} distinct positive triangle invariants.")
    print(f"PASS: {progression} controlled-Vandermonde progression cases.")
    print("PASS: exact cofactors, complement invariance, and radius inequality.")


if __name__ == "__main__":
    main()
