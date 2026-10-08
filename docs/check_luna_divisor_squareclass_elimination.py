#!/usr/bin/env python3
"""Exact checks for the same-squareclass divisor elimination identity."""
from math import isqrt


def divisors(h):
    return [n for n in range(1, h + 1) if h % n == 0]


def squareclass(n):
    d = n
    p = 2
    while p * p <= d:
        while d % (p * p) == 0:
            d //= p * p
        p += 1
    return d


def main():
    checked = 0
    for h in range(1, 81):
        for b0 in range(1, 41):
            for A in range(1, 401):
                sols = []
                for n in divisors(h):
                    v = 2 * A * n - n * n - b0 * b0
                    if v >= 0 and isqrt(v) ** 2 == v:
                        sols.append((n, isqrt(v)))
                by_class = {}
                for n, y in sols:
                    by_class.setdefault(squareclass(n), []).append((n, y))
                for d, rows in by_class.items():
                    for n, y in rows:
                        for np, yp in rows:
                            if n >= np:
                                continue
                            a2, b2 = n // d, np // d
                            a, bb = isqrt(a2), isqrt(b2)
                            assert d * a * a == n and d * bb * bb == np
                            lhs = (bb * y - a * yp) * (bb * y + a * yp)
                            rhs = (a * a - bb * bb) * (b0 * b0 - n * np)
                            assert lhs == rhs
                            checked += 1
                if A > (h + b0) ** 6:
                    for rows in by_class.values():
                        assert len(rows) <= 2
                        if len(rows) == 2:
                            assert rows[0][0] * rows[1][0] == b0 * b0
    print(f"PASS: {checked} exact pair identities and large-A squareclass checks")


if __name__ == "__main__":
    main()
