"""Exact checks supplementing the binary endpoint self-map obstruction."""

from fractions import Fraction as F
from itertools import product
from math import gcd, isqrt

from check_gaussian_reflection_replacement import gnorm
from check_mobius_actual_anchor import reduced_coefficients, minimal_matrix
from check_parabolic_shear_anchor_content import constants
from check_mobius_reciprocal_stretch_grid import realize


def multiply(A,B):
    a,b,c,d = A
    e,f,g,h = B
    return a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h


def finite_order_audit():
    count = 0
    seen = set()
    for M in product(range(-3,4),repeat=4):
        a,b,c,d = M
        det = a*d-b*c
        if det <= 0 or gcd(*M) != 1:
            continue
        P = (1,0,0,1)
        for n in range(1,13):
            P = multiply(P,M)
            if P[1] == P[2] == 0 and P[0] == P[3]:
                assert n in (1,2,3,4,6)
                if n > 1:
                    expected = {2:0,3:1,4:2,6:3}[n]
                    assert F((a+d)**2,det) == expected
                count += 1
                seen.add(n)
                break
    assert seen == {1,2,3,4,6}
    return count


def shear_core_audit():
    count = 0
    for a in range(1,31):
        for b in range(-60,61):
            if b == 0 or gcd(a,b) != 1:
                continue
            involution = (a,b,0,-a)
            assert multiply(involution,involution) == (a*a,0,0,a*a)
            shear = (a,b,0,a)
            assert multiply((1,0,0,-1),involution) == shear
            u,v = reduced_coefficients(shear)
            minimal,epsilon = minimal_matrix(u,v)
            Nu,Nv = gnorm(u),gnorm(v)
            assert gcd(Nu,Nv) == 1
            assert epsilon == (1 if b % 2 else 4)
            assert epsilon*abs(Nu-Nv) == 4*a*a
            assert abs(minimal[0]*minimal[3]-minimal[1]*minimal[2]) == a*a
            count += 1
    return count


def exponent_audit():
    count = 0
    for m in (*range(128,257),512,1024,4096):
        q,u,_,_ = constants(m)
        Fm = (m-1)**2//4
        h = F(m,8*Fm)
        eta = 1-q-F(3,2)*u-3*h
        assert q <= F(8,m)
        assert u <= F(1,2)+F(2,m)
        assert 3*h <= F(2,m)
        assert F(m*(m-1),Fm) <= 5
        assert eta >= F(1,4)-F(13,m) >= F(1,8)
        exponent = m*q+F(3,2)*u+F(3*m*(m-1),2*Fm)
        assert exponent <= F(33,2)
        # 40*2^(33/2)<2^22, certified after squaring.
        assert 40**2 * 2**33 < 2**44
        assert abs(float(m*(F(1,4)-eta))-7) < .1
        count += 1
    return count


def conformal_fixed_content_audit():
    count = 0
    for maximum in range(1,13):
        rows = [(1,0)] + [(1,s*k) for k in range(1,maximum+1) for s in (-1,1)]
        points,N = realize(rows)
        root = isqrt(N)
        assert root*root == N
        assert gcd(*points[0]) == root
        anchor = points[0]
        # The physical fixed point must be an axis point in a primitive realization.
        assert anchor[0]*anchor[1] == 0
        sign = 1 if anchor[1] == 0 else -1
        assert {(sign*x,-sign*y) for x,y in points} == set(points)
        count += 1
    for e in range(1,30):
        for e0 in range(e+1):
            symmetric = -e0 == -(e-e0)
            assert symmetric == (2*e0 == e)
            if symmetric:
                assert min(e0,e-e0)*2 == e
                count += 1
    return count


if __name__ == '__main__':
    orders = finite_order_audit()
    cores = shear_core_audit()
    exponents = exponent_audit()
    contents = conformal_fixed_content_audit()
    print(f'PASS: {orders} finite-order matrix cases; {cores} exact involution/core '
          f'normalizations; {exponents} exponent certificates; '
          f'{contents} conformal fixed-content cases.')
