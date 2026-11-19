"""Exact count of the primitive near-square divisors in Section 2 of answer.md.

Enumerates all odd parents with split prime factors up to the fixed limit.
The upper polylogarithmic cutoff is automatic here for B=4.
"""
from collections import Counter
from itertools import combinations
from math import gcd, isqrt


def split_only(k):
    if k % 2 == 0:
        return False
    q = 3
    while q * q <= k:
        if k % q == 0:
            if q % 4 != 1:
                return False
            while k % q == 0:
                k //= q
        q += 2
    return k % 4 == 1


def census(limit):
    counts = bytearray(limit + 1)
    # When N = dm, the angular cutoff becomes 256 t^4 m <= d.
    # In particular 256 m^2 <= dm <= limit.
    allowed = [m for m in range(1, isqrt(limit // 256) + 1) if split_only(m)]
    for s in range(1, isqrt(limit) + 1):
        bound = min(isqrt(limit - s*s), isqrt(2*s), isqrt(isqrt(limit // 256)))
        for t in range(1, bound + 1):
            if (s + t) % 2 == 0 or gcd(s, t) != 1:
                continue
            d = s*s + t*t
            largest_m = min(limit // d, d // (256 * t**4))
            for m in allowed:
                if m > largest_m:
                    break
                counts[d*m] += 1
    hist = Counter(counts)
    best = max(hist)
    print('Upper limit:', limit)
    print('Counts for nonempty parents:', sorted((n, v) for n, v in hist.items() if n))
    triangles = 0
    for parent, count in enumerate(counts):
        if count < 3:
            continue
        witnesses = []
        for m in allowed:
            if parent % m:
                continue
            d = parent // m
            s = isqrt(d)
            t = isqrt(d - s*s)
            if t and s*s + t*t == d and gcd(s, t) == 1 and 256*t**4*m <= d:
                witnesses.append((d, s, t))
        assert len(witnesses) == count
        for first, second, third in combinations(witnesses, 3):
            product = first[0]*second[0]*third[0]
            triangles += isqrt(product)**2 == product
        if count == best:
            print('Maximizing parent and its (d, s, t):', parent, witnesses)
    print('Unordered triples of divisors with square product:', triangles)


if __name__ == '__main__':
    census(10_000_000)
