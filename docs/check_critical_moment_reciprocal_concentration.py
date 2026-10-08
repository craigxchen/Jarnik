"""Exact reciprocal coefficient, concentration, and Gram checks."""

from fractions import Fraction as Q
from itertools import combinations


def det3(a):
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            - a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            + a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))


if __name__ == "__main__":
    # Normalize p2=1 and p1=1. Degree is at most one in m and two
    # in rho, so this grid verifies the entire determinant identity.
    for m in (0, 1):
        for rho in (Q(0), Q(1), Q(2)):
            p3 = (1+3*rho)/4
            gram = [[m, 1, 1], [1, 1, p3], [1, p3, rho]]
            assert det3(gram) == (1-rho)*(9*m*rho-m-8)/16

    roots = [Q(x) for x in (1, 5, 9, -7, -8)]
    assert sum(roots) == sum(x**3 for x in roots) == 0
    y = [1/x for x in roots]
    p = [Q(len(y))] + [sum(x**k for x in y) for k in range(1, 5)]
    for k in (2, 4):
        elementary = Q(0)
        for subset in combinations(y, k):
            value = Q(1)
            for x in subset:
                value *= x
            elementary += value
        assert elementary == 0
    assert p[1]**2 == p[2]
    assert 4*p[1]*p[3] == p[2]**2+3*p[4]
    rho = p[4]/p[2]**2
    gram = [[p[i+j] for j in range(3)] for i in range(3)]
    assert det3(gram) == p[2]**3*(1-rho)*(9*len(y)*rho-len(y)-8)/16 > 0
    assert rho > Q(len(y)+8, 9*len(y))
    assert p[2] < 4*max(x*x for x in y)
    print("Passed exact reciprocal coefficient identities, Gram determinant, and concentration fixture.")
