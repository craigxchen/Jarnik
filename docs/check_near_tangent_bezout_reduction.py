"""Exact coordinate and rectangle checks; the count bound has a tiling proof."""
from math import gcd


def complement(p, q):
    for r in range(-25, 26):
        for s in range(-25, 26):
            if p * s - q * r == 1:
                return r, s
    raise AssertionError((p, q))


count = 0
for p in range(-12, 13):
    for q in range(-12, 13):
        if gcd(p, q) != 1:
            continue
        r, s = complement(p, q)
        a, b = p*p + q*q, p*r + q*s
        assert a * (r*r + s*s) - b*b == 1
        assert (b*b + 1) % a == 0
        for xx in range(-3, 4):
            for yy in range(-3, 4):
                x, y = p*xx + q*yy, p*yy - q*xx
                assert (x - b*y) % a == 0
                assert p*x - q*y == a*xx
                assert q*x + p*y == a*yy
                assert x*x + y*y == a*(xx*xx + yy*yy)
                count += 1
        for m in range(-3, 4):
            for y in range(-3, 4):
                x = a*m + b*y
                assert (p*x - q*y) % a == 0
                assert (q*x + p*y) % a == 0
print('PASS:', count, 'exact Bezout coordinates and norm comparisons')

for p in range(2, 301):
    a, h = p*p + 1, p//2
    k, count = p*h, 0
    for y in range(h+1):
        lower = -((k-p*y)//a)  # ceil((-k+p*y)/a)
        upper = (p*y)//a
        assert lower == upper == 0
        count += upper-lower+1
    assert count == h+1
    assert 8*(a*count-k*h) >= p*a
print('PASS: 299 exact rectangle counts and linear area discrepancy')
