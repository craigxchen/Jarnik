"""Exact covariance and common-content checks for the discriminant lift."""
from fractions import Fraction as F
from itertools import product
from math import gcd, lcm

from check_gaussian_reflection_replacement import gmul, gconj, gnorm
from check_mobius_conductor_transfer import fixtures


def mmul(M, L):
    a, b, c, d = M
    e, f, g, h = L
    return (a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h)


def scalar_matrix(z):
    x, y = z
    return (x, -y, y, x)


def primitive(row):
    common = gcd(*row)
    return tuple(x//common for x in row)


def invariants(M):
    a, b, c, d = M
    T = a*a+b*b+c*c+d*d
    delta = a*d-b*c
    Gamma = (b*b+d*d-a*a-c*c, 2*(a*b+c*d))
    assert gnorm(Gamma) == T*T-4*delta*delta
    return T, delta, Gamma


def canonical(M, rows):
    T, delta, Gamma = invariants(M)
    points = []
    h = 1
    for H in rows:
        z = gmul(Gamma, gmul(H, H))
        f = gnorm(H)
        point = tuple(F(x, f) for x in z)
        points.append(point)
        for x in point:
            h = lcm(h, x.denominator)
    lifted = [tuple(int(h*x) for x in z) for z in points]
    content = gcd(*(x for z in lifted for x in z))
    assert gcd(content, h) == 1
    r = F(h, content)
    canon = [tuple(x//content for x in z) for z in lifted]
    center = r*T
    defect_root = 2*r*delta
    for H, z in zip(rows, canon):
        assert center*center-gnorm(z) == defect_root*defect_root
        a, b, c, d = M
        x, y = H
        stretch = F((a*x+b*y)**2+(c*x+d*y)**2, gnorm(H))
        assert center-z[0] == 2*r*stretch
    return canon, center, defect_root, h, content


def main():
    sources = fixtures()
    sources.append([(1,0),(1,1),(1,-1)])
    covariance_count = content_count = nontrivial = 0
    matrices = []
    for M in product(range(-1,2), repeat=4):
        T, delta, Gamma = invariants(M)
        if delta and gnorm(Gamma):
            matrices.append(M)
    matrices.extend([(1,155,0,1),(5,1,0,2),(3,4,-2,1)])
    for M, rows in product(matrices, sources):
        base = canonical(M, rows)
        T, delta, Gamma = invariants(M)
        for Hj in rows:
            reanchored = [primitive(gmul(H, gconj(Hj))) for H in rows]
            for lam in [(1,0),(1,1),(2,-1)]:
                moved = mmul(scalar_matrix(lam), mmul(M,scalar_matrix(Hj)))
                scale = gnorm(lam)*gnorm(Hj)
                newT, newdelta, newGamma = invariants(moved)
                assert newT == scale*T and newdelta == scale*delta
                expected = gmul(Gamma,gmul(Hj,Hj))
                assert newGamma == tuple(gnorm(lam)*x for x in expected)
                result = canonical(moved,reanchored)
                assert result[:3] == base[:3]
                common = gcd(*moved)
                reduced = tuple(x//common for x in moved)
                assert canonical(reduced,reanchored)[:3] == base[:3]
                covariance_count += 1
        # Orthogonal reflection of the output preserves Gamma and center.
        a,b,c,d = M
        reflected = (a,b,-c,-d)
        result = canonical(reflected,rows)
        assert result[0:2] == base[0:2]
        assert result[2] == -base[2]

    for A,B,D in product(range(1,7),range(-6,7),range(-6,7)):
        if not D or gcd(A,B,D) != 1:
            continue
        M = (A,B,0,D)
        T,delta,Gamma = invariants(M)
        if not gnorm(Gamma):
            continue
        for rows in sources:
            canon,center,defect_root,h,content = canonical(M,rows)
            d0 = gcd(content,h*T)
            assert d0 == gcd(content,T)
            assert 2*A % d0 == 0
            # The smallest integer dilation of the canonical tuple whose
            # center is integral is exactly the center denominator.
            multiplier = content//d0
            assert center.denominator == multiplier
            assert all((multiplier*x).denominator == 1 for x in (center,defect_root))
            content_count += 1
            nontrivial += d0 > 1
    print(f'PASS: {covariance_count} exact actual-anchor/output-factor covariance checks.')
    print('PASS: canonical integer tuple, rational center, signed defect root, and stretch/radial identity.')
    print(f'PASS: {content_count} primitive triangular cases; integer-center content divides 2A ({nontrivial} nontrivial).')


if __name__ == '__main__':
    main()
