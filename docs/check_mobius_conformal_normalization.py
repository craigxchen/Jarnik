"""Exact finite checks of output-conformal minimum and padding invariance."""
from itertools import product
from math import gcd
from check_gaussian_reflection_replacement import ggcd, gdivexact, gmul, gnorm
from check_mobius_conductor_transfer import fixtures, least_norm


def matrix_from_doubled(alpha2, beta2):
    x, y = alpha2
    u, v = beta2
    entries2 = (x+u, -y+v, y+v, x-u)
    if any(e % 2 for e in entries2):
        return None
    return tuple(e//2 for e in entries2)


def invariants(matrix):
    a, b, c, d = matrix
    T = a*a+b*b+c*c+d*d
    delta = a*d-b*c
    return T, T*T-4*delta*delta


def good_determinant(matrix):
    a,b,c,d = matrix
    determinant = abs(a*d-b*c)
    _,K = invariants(matrix)
    while gcd(determinant,K) > 1:
        determinant //= gcd(determinant,K)
    return determinant


def image(rows, matrix):
    a, b, c, d = matrix
    out = []
    for x, y in rows:
        u, v = a*x+b*y, c*x+d*y
        g = gcd(u, v)
        out.append((u//g, v//g))
    return out


def main():
    coefficients = [(x,y) for x,y in product(range(-3,4), repeat=2) if x or y]
    multipliers = [(x,y) for x,y in product(range(-3,4), repeat=2) if x or y]
    cases = 0
    for u, v in product(coefficients, repeat=2):
        if gnorm(ggcd(u,v)) != 1 or gnorm(u) == gnorm(v):
            continue
        D = ggcd((2,0), (u[0]-v[0],u[1]-v[1]))
        gamma = gdivexact((2,0), D)
        epsilon = 4//gnorm(D)
        assert gnorm(gamma) == epsilon and epsilon in (1,2,4)
        matrix = matrix_from_doubled(gmul(gamma,u), gmul(gamma,v))
        assert matrix is not None and gcd(*matrix) == 1
        T, K = invariants(matrix)
        best_determinant = good_determinant(matrix)
        assert 2*T == epsilon*(gnorm(u)+gnorm(v))
        assert K == epsilon**2*gnorm(u)*gnorm(v)
        for eta in multipliers:
            candidate = matrix_from_doubled(gmul(eta,u),gmul(eta,v))
            divisible = (gmul(eta,(u[0]-v[0],u[1]-v[1]))[0] % 2 == 0
                         and gmul(eta,(u[0]-v[0],u[1]-v[1]))[1] % 2 == 0)
            assert (candidate is not None) == divisible
            if candidate is not None:
                t,k = invariants(candidate)
                assert gnorm(eta) % epsilon == 0
                assert t >= T and k >= K
                assert t*epsilon == T*gnorm(eta)
                assert k*epsilon**2 == K*gnorm(eta)**2
                assert best_determinant % good_determinant(candidate) == 0
        cases += 1
    padding = 0
    for rows in fixtures():
        base = least_norm(image(rows,(2,0,0,1)))
        for n in range(1,31):
            matrix = (2*n,-1,2,n)
            T,K = invariants(matrix)
            assert gcd(*matrix) == 1
            assert T == 5*(n*n+1) and K == 9*(n*n+1)**2
            assert least_norm(image(rows,matrix)) == base
            padding += 1
    print(f'PASS: {cases} coprime coefficient pairs, exact integrality and simultaneous minimum.')
    print(f'PASS: {padding} primitive conformal-padding cases preserve least radius.')
    print('PASS: minimal output representative also maximizes the good determinant part.')


if __name__ == '__main__':
    main()
