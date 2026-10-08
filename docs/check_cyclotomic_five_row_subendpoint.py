#!/usr/bin/env python3
"""Exact checks for cyclotomic_five_row_subendpoint_family.md.

Dependency-free. The default k=19, d=1,3 verifies the smallest
coprime-order fixture. Use --d 1 3 5 for an extended check and
--k 21 for the fixture in the note.
All circle and arc checks use integer arithmetic; no floating-point
angle calculation enters the assertions.
"""

from __future__ import annotations

import argparse
from itertools import combinations
from math import gcd


def factor(n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def divisors(n: int) -> list[int]:
    ds = [1]
    for p, a in factor(n).items():
        old = list(ds)
        q = 1
        for _ in range(a):
            q *= p
            ds.extend(q * d for d in old)
    return sorted(ds)


def phi(n: int) -> int:
    ans = n
    for p in factor(n):
        ans -= ans // p
    return ans


def mobius(n: int) -> int:
    if any(a > 1 for a in factor(n).values()):
        return 0
    return -1 if len(factor(n)) % 2 else 1


def fib_pair(n: int) -> tuple[int, int]:
    """Return (F_n,F_{n+1}) by exact fast doubling."""
    if n == 0:
        return 0, 1
    a, b = fib_pair(n // 2)
    c = a * (2 * b - a)
    d = a * a + b * b
    return (d, c + d) if n & 1 else (c, d)


G = tuple[int, int]


def mul(a: G, b: G) -> G:
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def conj(a: G) -> G:
    return a[0], -a[1]


def sub(a: G, b: G) -> G:
    return a[0] - b[0], a[1] - b[1]


def norm(a: G) -> int:
    return a[0] * a[0] + a[1] * a[1]


def nearest(n: int, d: int) -> int:
    return (2 * n + d) // (2 * d)


def divmod_gaussian(a: G, b: G) -> tuple[G, G]:
    den = norm(b)
    q = (
        nearest(a[0] * b[0] + a[1] * b[1], den),
        nearest(a[1] * b[0] - a[0] * b[1], den),
    )
    return q, sub(a, mul(q, b))


def gcd_gaussian(a: G, b: G) -> G:
    while b != (0, 0):
        _, r = divmod_gaussian(a, b)
        assert norm(r) < norm(b)
        a, b = b, r
    return a


def exact_div(a: G, b: G) -> G:
    q, r = divmod_gaussian(a, b)
    assert r == (0, 0)
    return q


def g_fib(odd_index: int) -> G:
    assert odd_index > 0 and odd_index & 1
    r = (odd_index - 1) // 2
    f_r, f_r1 = fib_pair(r)
    return f_r1, f_r


def orders(k: int) -> list[int]:
    x, y, z, w = k - 2, k, k + 2, k + 4
    assert k >= 19 and k & 1 and k % 3 in (0, 1)
    assert all(gcd(a, b) == 1 for a, b in combinations((x, y, z, w), 2))
    return [x * z * w, y * y * w, y * z * z, y * z * w]


def check_layers(k: int) -> tuple[list[int], int, int]:
    ns = orders(k)
    h = ns[0]
    assert h == min(ns) and max(ns) < 3 * h
    E = sorted(set().union(*(set(divisors(n)) for n in ns)))
    assert [e for e in E if e >= h] == sorted(ns)
    U = sum(phi(e) for e in E)
    assert U == 4 * k**3 + 15 * k**2 - 4 * k - 24

    rows = [
        [int(n % e == 0) for e in E] for n in ns
    ]
    rows.append([sum(rows[j][i] for j in (0, 1, 2)) for i in range(len(E))])
    base_index = E.index(1)
    assert [row[base_index] for row in rows] == [1, 1, 1, 1, 3]

    L = sum(phi(e) * (max(row[i] for row in rows) -
                       min(row[i] for row in rows))
            for i, e in enumerate(E))
    assert L == U + 3 * k + 4
    assert 4 * h - L == k * k - 15 * k - 44 > 0

    # Möbius transform of each regular-polygon incidence vector.
    for j, row in enumerate(rows):
        for d in E:
            val = sum(mobius(e // d) * row[i]
                      for i, e in enumerate(E) if e % d == 0)
            expected = int(d == ns[j]) if j < 4 else sum(
                int(d == ns[t]) for t in (0, 1, 2)
            )
            assert val == expected, (j, d, val, expected)

    # Independent Ramanujan-sum moment audit for every odd t through h.
    # The sum over e|n of c_e(t) equals n if n|t, and zero otherwise.
    terms = [
        [(q, q * mobius(e // q)) for q in divisors(e)]
        for e in E
    ]
    for t in range(1, h + 1, 2):
        c = [sum(weight for q, weight in qs if t % q == 0)
             for qs in terms]
        for j, row in enumerate(rows):
            got = sum(a * ce for a, ce in zip(row, c))
            expected = (
                sum(ns[s] for s in (0, 1, 2) if t % ns[s] == 0)
                if j == 4 else (ns[j] if t % ns[j] == 0 else 0)
            )
            assert got == expected, (j, t, got, expected)
    print(f"k={k}: h={h}, |E|={len(E)}, U={U}, L={L}, margin={4*h-L}")
    return ns, h, L


def check_gaussian(ns: list[int], d: int) -> None:
    assert d > 0 and d & 1
    H = [g_fib(n * d) for n in ns]
    F, Fb = (1, 2), (1, -2)
    patterns = ((0,), (1,), (2,), (3,), (0, 1, 2))
    Z: list[G] = []
    for j, forward in enumerate(patterns):
        val = Fb if j == 4 else F
        for t in range(4):
            val = mul(val, H[t] if t in forward else conj(H[t]))
        Z.append(val)
    N = norm(Z[0])
    assert all(norm(z) == N for z in Z)
    assert len(set(Z)) == 5

    D = Z[0]
    for z in Z[1:]:
        D = gcd_gaussian(D, z)
    P = [exact_div(z, D) for z in Z]
    primitive_norm = norm(P[0])
    assert all(norm(p) == primitive_norm for p in P)
    assert len(set(P)) == 5
    # Check the gcd of all primitive rows, not just a chosen triple.
    check_gcd = P[0]
    for p in P[1:]:
        check_gcd = gcd_gaussian(check_gcd, p)
    assert norm(check_gcd) == 1

    # Pairwise positive dot products put all arguments in one minor arc.
    assert all((a[0] * b[0] + a[1] * b[1]) > 0
               for a, b in combinations(P, 2))
    chord2 = max(norm(sub(a, b)) for a, b in combinations(P, 2))
    # If theta<=pi, arc/chord<=pi/2. The rational inequality
    # 7*chord^4<R^2 implies arc/sqrt(R)<pi/(2*7^(1/4))<1.
    assert 7 * chord2 * chord2 < primitive_norm
    print(
        f"d={d}: primitive norm bits={primitive_norm.bit_length()}, "
        f"gcd norm bits={norm(D).bit_length()}, "
        f"7*chord^4/R^2<1"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=19)
    parser.add_argument("--d", type=int, nargs="*", default=[1, 3])
    args = parser.parse_args()
    ns, _, _ = check_layers(args.k)
    for d in args.d:
        check_gaussian(ns, d)


if __name__ == "__main__":
    main()
