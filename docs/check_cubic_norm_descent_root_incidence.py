"""Exact checks for the cubic norm-descent construction and incidence bound."""

from fractions import Fraction as F
from itertools import combinations


def trim(a):
    a = list(map(F, a))
    while a and not a[-1]:
        a.pop()
    return tuple(a)


def add(a, b):
    return trim([(a[j] if j < len(a) else 0) +
                 (b[j] if j < len(b) else 0)
                 for j in range(max(len(a), len(b)))])


def scale(a, c):
    return trim([c*x for x in a])


def mul(a, b):
    if not a or not b:
        return ()
    out = [F(0)] * (len(a)+len(b)-1)
    for j, x in enumerate(a):
        for k, y in enumerate(b):
            out[j+k] += x*y
    return trim(out)


def at_i(a):
    re = sum(((-1)**(j//2))*x for j, x in enumerate(a) if j % 2 == 0)
    im = sum(((-1)**(j//2))*x for j, x in enumerate(a) if j % 2 == 1)
    return re, im


def check_quartic(delta, r):
    u = 2*r/(1+delta*r*r)
    v0 = (1-delta*r*r)/(1+delta*r*r)
    assert u and v0 and v0*v0+delta*u*u == 1
    k = v0/(4*delta*u*u)
    s = (F(1), F(0), F(1))
    a = add((0, F(-3, 2), 0, F(-1, 2)), scale(mul(s, s), k))
    b = mul(s, (k*(1-2*delta*u*u), F(-1, 2), k))
    c = mul(s, (-(delta*u*u+1)/(2*delta*u), 2*k*u))
    assert add(mul(a, a), (1,)) == add(mul(b, b), scale(mul(c, c), delta))
    assert at_i(a) == (0, -1)
    assert at_i(b) == at_i(c) == (0, 0)
    assert len(a) == len(b) == 5 and a[-1] == b[-1] == k
    assert len(c) == 4


def main():
    deltas = tuple(map(F, (-7, -3, -2, -1))) + (F(-3, 2), F(-2, 3),
              F(-1, 2), F(-1, 3), F(1, 3), F(1, 2), F(2, 3),
              F(1), F(3, 2), F(2), F(3), F(7))
    for delta in deltas:
        check_quartic(delta, F(1, 5))

    # Three conjugate values p_j satisfy e2(p)=2i sum(p).
    # If p_1=p_2=0, then p_3=0; no numerical fixture is needed.
    disjoint_counts = {}
    for size in (2, 3):
        disjoint_counts[size] = sum(
            all(a & b == 0 for a, b in combinations(cuts, 2))
            for cuts in combinations(range(1, 16), size)
        )
    assert disjoint_counts == {2: 25, 3: 10}

    maximum = 0
    for bits in range(1 << 14):
        masks = [m for m in range(1, 15) if bits >> (m-1) & 1] + [15]
        if all(sum(bool(m & (1 << j)) for m in masks) <= 3 for j in range(4)):
            maximum = max(maximum, len(masks))
    assert maximum == 7
    print('PASS: 16 exact quartic identities, fixed-root values, and leading coefficients.')
    print('PASS: disjoint-cut orbit patterns and at most seven Gaussian-rational labels.')


if __name__ == '__main__':
    main()
