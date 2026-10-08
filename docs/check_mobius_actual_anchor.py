"""Exact tests for normalization at actual source anchors."""
from itertools import product
from math import gcd

from check_gaussian_reflection_replacement import gmul, gconj, gnorm, ggcd, gdivexact
from check_mobius_conductor_transfer import fixtures, least_norm
from check_mobius_conformal_normalization import matrix_from_doubled, invariants, good_determinant, image


def primitive(z):
    r = gcd(*z)
    return (z[0]//r,z[1]//r)


def reduced_coefficients(matrix):
    a,b,c,d = matrix
    alpha2 = (a+d,c-b)
    beta2 = (a-d,c+b)
    common = ggcd(alpha2,beta2)
    return gdivexact(alpha2,common),gdivexact(beta2,common)


def minimal_matrix(u,v):
    common = ggcd((2,0),(u[0]-v[0],u[1]-v[1]))
    gamma = gdivexact((2,0),common)
    epsilon = gnorm(gamma)
    matrix = matrix_from_doubled(gmul(gamma,u),gmul(gamma,v))
    assert matrix is not None and gcd(*matrix) == 1
    return matrix,epsilon


def multiply_input(matrix,h):
    a,b,c,d = matrix
    x,y = h
    return (a*x+b*y,-a*y+b*x,c*x+d*y,-c*y+d*x)


def good_part(n,k):
    out = abs(n)
    while gcd(out,k) > 1:
        out //= gcd(out,k)
    return out


def main():
    cases = 0
    scaled_cases = 0
    for entries in product(range(-2,3),repeat=4):
        a,b,c,d = entries
        if not a*d-b*c or gcd(*entries) != 1:
            continue
        T,K = invariants(entries)
        if K == 0:
            continue
        u,v = reduced_coefficients(entries)
        minimal,epsilon = minimal_matrix(u,v)
        U,V = gnorm(u),gnorm(v)
        k0 = gcd(U,V)
        aa,bb = V//k0,U//k0
        assert k0 % 2 == 1 and gcd(aa,bb) == 1
        T0,K0 = invariants(minimal)
        delta = minimal[0]*minimal[3]-minimal[1]*minimal[2]
        assert T0 % k0 == 0 and K0 % (k0*k0) == 0 and delta % k0 == 0
        Tbase,Kbase,deltabase = T0//k0,K0//(k0*k0),delta//k0
        Dbase = good_part(deltabase,Kbase)
        for rows in fixtures():
            source_norm = least_norm(rows)
            target_norm = least_norm(image(rows,minimal))
            for j,h in enumerate(rows):
                new_rows = [primitive(gmul(z,gconj(h))) for z in rows]
                assert new_rows[j] == (1,0)
                assert all(gnorm(z) % 2 and gcd(*z) == 1 for z in new_rows)
                assert least_norm(new_rows) == source_norm
                pre_u,pre_v = gmul(u,h),gmul(v,gconj(h))
                common = ggcd(pre_u,pre_v)
                uj,vj = gdivexact(pre_u,common),gdivexact(pre_v,common)
                factors = gmul(ggcd(u,gconj(h)),ggcd(v,h))
                assert gnorm(factors) == gnorm(common)
                Uj,Vj = gnorm(uj),gnorm(vj)
                assert Uj*gnorm(common) == U*gnorm(h)
                assert Vj*gnorm(common) == V*gnorm(h)
                kj = gcd(Uj,Vj)
                assert Uj == bb*kj and Vj == aa*kj
                moved,epsj = minimal_matrix(uj,vj)
                assert epsj == epsilon
                tj,kkj = invariants(moved)
                dj = moved[0]*moved[3]-moved[1]*moved[2]
                assert (tj,kkj,dj) == (Tbase*kj,Kbase*kj*kj,deltabase*kj)
                assert good_determinant(moved) == good_part(Dbase,kj)
                raw_moved = multiply_input(minimal,h)
                assert least_norm(image(new_rows,raw_moved)) == target_norm
                assert least_norm(image(new_rows,moved)) == target_norm
                cases += 1
        scaled_cases += 1
    print(f'PASS: {cases} actual-anchor cases across {scaled_cases} primitive nonconformal matrices.')
    print('PASS: source and full target least norms, exact Gaussian contents, constant epsilon, and all T/K/determinant formulas.')


if __name__ == '__main__':
    main()
