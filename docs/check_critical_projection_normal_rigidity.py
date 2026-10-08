"""Exact normal covariance and algebra; no synthetic endpoint tuples asserted."""

from fractions import Fraction as F
from itertools import product
from math import gcd

from check_gaussian_reflection_replacement import gmul, gconj, gnorm
from check_mobius_conductor_transfer import fixtures
from check_mobius_reciprocal_stretch_grid import realize
from check_parabolic_shear_anchor_content import constants


def multiply(A,B):
    a,b,c,d = A
    e,f,g,h = B
    return a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h


def gram(M):
    a,b,c,d = M
    return a*a+c*c-b*b-d*d,2*(a*b+c*d)


def normal(C):
    g = gcd(*C)
    assert g > 0
    return tuple(x//g for x in C)


def covariance_audit():
    cases = 0
    for rows in fixtures():
        points,N = realize(rows)
        for M in product(range(-2,3),repeat=4):
            a,b,c,d = M
            if a*d-b*c == 0 or gram(M) == (0,0):
                continue
            W = gram(M)
            C = gmul(W,points[0])
            for H,z in zip(rows,points):
                x,y = H
                R = (x,-y,y,x)
                moved = multiply(M,R)
                u,v = moved[0],moved[2]
                L = (u,v,-v,u)
                raw = multiply(L,moved)
                content = gcd(*raw)
                triangular = tuple(q//content for q in raw)
                assert triangular[2] == 0
                Cnew = gmul(gram(triangular),z)
                scale = F((u*u+v*v)*gnorm(H),content*content)
                assert Cnew == tuple(scale*q for q in C)
                assert scale > 0 and normal(Cnew) == normal(C)
                cases += 1
            for lam in ((1,1),(2,-1),(-1,3)):
                x,y = lam
                left = (x,-y,y,x)
                padded = multiply(left,M)
                assert gram(padded) == tuple(gnorm(lam)*q for q in W)
                reflected = multiply((1,0,0,-1),padded)
                assert gram(reflected) == gram(padded)
    return cases


def inverse_and_product_audit():
    inverse_cases = product_cases = 0
    for M in product(range(-4,5),repeat=4):
        a,b,c,d = M
        if a*d-b*c == 0:
            continue
        W = gram(M)
        if W == (0,0):
            continue
        J = (d,-b,-c,a)
        Wj = gram(J)
        same_oriented_axis = normal(W) == normal(Wj)
        assert same_oriented_axis == (a+d == 0)
        if same_oriented_axis:
            square = multiply(M,M)
            assert square[1] == square[2] == 0 and square[0] == square[3]
            inverse_cases += 1
    matrices = [(0,b,c,0) for b,c in product(range(-8,9),repeat=2)
                if b*c and abs(c)>abs(b)]
    for A,B in product(matrices,repeat=2):
        assert normal(gram(A)) == normal(gram(B)) == (1,0)
        C = multiply(A,B)
        assert C[1] == C[2] == 0
        P = (1,0,0,1)
        finite_order = False
        for _ in range(6):
            P = multiply(P,C)
            if P[0] == P[3]:
                finite_order = True
                break
        assert finite_order == (abs(C[0]) == abs(C[3]))
        if finite_order:
            assert gram(C) == (0,0)
        product_cases += 1
    # Same axes alone are insufficient without the finite-order premise.
    assert gram(multiply((0,1,2,0),(0,1,3,0))) != (0,0)
    return inverse_cases,product_cases


def exponent_audit():
    count = 0
    for m in (*range(512,641),768,1024,2048,4096):
        q,u,h,_ = constants(m)
        Fm = (m-1)**2//4
        g = F(m,4*Fm)
        z = (1-q)/2-g
        b = h-g
        log2_D = 2+F(m*(m-1),Fm)+m*b
        nu = 2*h+u-F(1,2)+q
        mu = (h-z)/2
        assert log2_D <= 14
        assert h <= F(7,m) and q <= F(8,m) and m*q <= 8
        assert u <= F(1,2)+F(2,m) <= F(2,3)
        assert F(1,2)-F(6,m) <= z <= F(1,2)
        assert z >= F(1,4)
        assert nu <= F(24,m)
        assert -F(1,4) <= mu <= -F(1,4)+F(13,2*m)
        assert 2*nu+mu <= -F(1,4)+F(109,2*m) <= -F(1,8)
        assert (2-F(4,m))*z-u >= F(1,4)
        assert 20+u+4*z < 23
        assert abs(float(m*(2*nu+mu+F(1,4)))-209/8) < .3
        count += 1
    # Height/angular constants and the stated uniform radius threshold.
    assert 5**3 * 2**110 < 2**117
    assert 128*2**7+2 < 2**15
    assert 2*39+15 == 93 and 8*93 == 744
    assert 744//4 > 23
    return count


def common_axis_composite_audit():
    count = 0
    pairs = [(r,s) for r in range(2,7) for s in range(1,r)]
    for r,s in pairs:
        for t,u in pairs:
            for x,y in ((1,0),(2,1),(-1,3)):
                R = (x,-y,y,x)
                LA = (2,-3,3,2)
                LB = (-1,-4,4,-1)
                for sign in (-1,1):
                    A = multiply(LA,multiply((r,0,0,s),R))
                    B = multiply((1,0,0,sign),multiply(LB,multiply((t,0,0,u),R)))
                    assert normal(gram(A)) == normal(gram(B))
                    a,b,c,d = A
                    C = multiply(B,(d,-b,-c,a))  # Scalar multiple of B A^{-1}.
                    aa,bb,cc,dd = C
                    trace = aa*aa+bb*bb+cc*cc+dd*dd
                    determinant = abs(aa*dd-bb*cc)
                    ratio = F(t*s,u*r)
                    kappa = max(ratio,1/ratio)
                    assert F(trace,determinant) == kappa+1/kappa
                    count += 1
    return count


if __name__ == '__main__':
    covariance = covariance_audit()
    inverse,products = inverse_and_product_audit()
    exponents = exponent_audit()
    composites = common_axis_composite_audit()
    print(f'PASS: {covariance} exact physical reanchoring covariances; '
          f'{inverse} inverse/common-oriented-axis cases; '
          f'{products} finite-order product cases; {composites} common-axis '
          f'composite identities; {exponents} exponent certificates.')
