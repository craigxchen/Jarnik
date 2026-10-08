"""Finite exact checks supporting the analytic uniform-tail proof."""

from fractions import Fraction as Q


def polynomial(roots):
    coefficients = [Q(1)]
    for z in roots:
        updated = [Q(0)]*(len(coefficients)+1)
        for i, c in enumerate(coefficients):
            updated[i] += c
            updated[i+1] -= c/z
        coefficients = updated
    return coefficients


def check_fixture(roots):
    roots = list(map(Q, roots))
    minimum = min(abs(z) for z in roots)
    roots = sorted((z/minimum for z in roots), key=abs)
    assert len({abs(z) for z in roots}) == len(roots)
    s = (len(roots)-1)//2
    assert all(sum(z**(2*r+1) for z in roots) == 0 for r in range(s))
    coefficients = polynomial(roots)
    assert coefficients[0] == 1
    assert all(coefficients[j] == 0 for j in range(2, len(coefficients), 2))
    y = [1/z for z in roots]
    p1, p2 = sum(y), sum(t*t for t in y)
    assert p1*p1 == p2 < 4
    assert all((j+1)*t*t <= 4 for j, t in enumerate(y))
    for K in range(len(y)+1):
        assert sum(abs(t)**3 for t in y[K:]) <= (
            max((abs(t) for t in y[K:]), default=Q(0))*sum(t*t for t in y[K:]))


if __name__ == "__main__":
    check_fixture([-3, 1, 2])
    check_fixture([1, 5, 9, -7, -8])
    check_fixture([-2, -10, -18, 14, 16])
    # The logarithmic remainder uses 1/k<=1/3 for every k>=3,
    # and the geometric tail bound 1/[3(1-r)]<=2/3 for 0<=r<=1/2.
    assert Q(1, 3)/(1-Q(1, 2)) == Q(2, 3)
    assert Q(2, 3)*2*4 == Q(16, 3)
    assert min(r*(4-t)+(4-r)*t for r in range(5) for t in range(5)
               if (r,t) not in ((0,0),(4,4))) == 4
    # For u>1, the last real-axis comparison has a positive square gap.
    for u in (Q(3, 2), Q(2), Q(5), Q(100)):
        assert u*u-(2*u-2) == (u-1)**2+1 > 0
    print("Passed finite critical products, reciprocal bounds, and analytic-estimate constants.")
